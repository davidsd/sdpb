# bigint_syrk_blas() algorithm description

## Introduction

Suppose that we want to calculate matrix square Q = P^T P, where P is a `BigFloat` matrix.

(In SDPB, P is `schur_off_diagonal` and Q is `Q`,
see [compute_Q.cxx](../step/initialize_schur_complement_solver/compute_Q.cxx))

How can we do that?

For double-precision arithmetic, the fastest way is BLAS routine `cblas_dsyrk()`.

For multiprecision (`BigFloat`), one can `El::Syrk()` from Elemental. But it doesn't benefit from heavy BLAS
optimizations and thus is much slower.

Since we are aiming for the best performance, let's find a way to use the power of highly optimized BLAS routines. How?

General scheme:

1. Split P into several double-precision matrices P_i
2. Compute Q_i := P_i^T P_i
3. Restore Q from Q_i

How can we do that? Let's see.

## Multiplying big integers: Chinese Remainder Theorem

If we have a multiprecision integer matrix, then the following trick works:

1. Choose a set of primes, such that each prime `p` satisfies a condition `p^2 * k < 2^53`, where `k = P.Height()`, and
   product of all primes exceeds `max(Q)`.
2. Calculate residues of each element modulo some primes. Since all primes are small enough, all calculations will
   fit into double precision (53 bits).
3. Convert each residue to double and calculate residue matrix Q_i := P_i^T P_i using `cblas_dsyrk`.
4. Restore Q from Q_i using [Chinese Remainder Theorem](https://en.wikipedia.org/wiki/Chinese_remainder_theorem) (CRT).

This algorithm is implemented in [FLINT library](https://flintlib.org/),
see [mul_blas.c](https://github.com/flintlib/flint2/blob/trunk/src/fmpz_mat/mul_blas.c). We reuse its parts in our code.

## From BigFloat to BigInt and back: matrix normalization

Can we utilize the same trick for BigFloats, without precision loss?

Yes! If we normalize columns of P and multiply P by 2^N, where N is a number of bits, do all calculations and then
restore Q.

So the algorithm is:

1. Normalize columns of P, multiply by 2^N. Let's call it P'.
2. Choose a set of primes.
3. Calculate P'_i - residues of P' modulo each prime, convert to double.
4. Calculate Q'_i := P'_i^T P'_i using `cblas_dsyrk`.
5. Restore Q' from residues Q'_i using CRT. Diagonal elements of Q' are all equal to 2^2N, thanks to P' normalization.
6. Restore Q from Q'.

Matrix types:

- P and Q are `Matrix<BigFloat>`
- P' and Q' are also `Matrix<BigFloat>`, in fact containing big integers (i.e. each element can be cast to
  FLINT `fmpz_t` type without losing precision).
- P'_i and Q'_i are `Matrix<double>`, in fact containing integers (can be cast to `uint32_t`).

Class [Matrix_Normalizer](Matrix_Normalizer.hxx) performs steps 1 and 6.

Class [Fmpz_Comb](Fmpz_Comb.hxx), which is a wrapper of FLINT's `fmpz_comb_t`,
helps us with steps 2, 3, and 5.

## Going parallel

In SDPB, everything is more complicated than described above, because P is not a `Matrix<BigFloat>`, but
a `vector<DistMatrix<BigFloat>>`.

This means that P is split into horizontal bands (or blocks), and each block is stored on one or several cores.

In turn, Q is `DistMatrix<BigFloat>` distributed over all cores.

In the old SDPB algorithm, each core (or group of cores, if `procGranularity > 1`) was calculating a contribution to Q (
called `Q_group`) from its blocks. Adding up all Q_groups, one got the resulting Q
(see [synchronize_Q.cxx](../step/initialize_schur_complement_solver/synchronize_Q.cxx)).

This is inefficient in terms of memory, since we allocate memory for `Q_group` for each group of cores.

Now we can avoid this duplication, if all cores on a single node are working together on the blocks,
using [MPI shared memory window](https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPI_Win_allocate_shared.3.html).

The parallel algorithm goes as follows (see [compute_Q](../step/initialize_schur_complement_solver/compute_Q.cxx)
and [bigint_syrk_blas](BigInt_Shared_Memory_Syrk_Context.cxx) functions):

1. Choose a set of primes. This is done only once at SDPB start,
   see [initialize_bigint_syrk_context](initialize_bigint_syrk_context.hxx).
2. Calculate column norms of P (requires synchronization over all P blocks).
3. Normalize each P block (divide by column norms) and multiply it by 2^N, where `N = El::gmp::Precision()`.
4. For each block on a node, calculate its residues and put them into shared memory window.
   The residues are stored in [Block_Residue_Matrices_Window](Block_Residue_Matrices_Window.hxx) in the following order:

```
block_1 mod prime_1
block_2 mod prime_1
block_3 mod prime_1
...
block_1 mod prime_2
block_2 mod prime_2
...
```

5. Each core takes part of this window and calls `cblas_dsyrk()`.
   Currently, we use a simple round-robin algorithm:
   core 1 takes a matrix composed of all blocks mod prime_1, core 2 - all blocks mode prime_2, etc.
   (TODO: if we have more cores than primes, then we should split a job for each prime into several cores).
   The results are stored in [Residue_Matrices_Window](Residue_Matrices_Window.hxx):

```
Q_group mod prime_1
Q_group mod prime_2
Q_group mod prime_3
...
```

6. `Q_group` is now a `DistMatrix<BigFloat>` distributed over all cores of a node. If some element `Q_group(i,j)` is
   stored e.g. on core 1, this core restores its value from the `Residue_Matrices_Window` using CRT.
7. Restore Q_group (see `Matrix_Normalizer.restore_Q`), i.e. divide by 2^2N and remove normalization.
8. Calculate global Q, which is `DistMatrix<BigFloat>` distributed over all cores, as a sum of all Q_groups.
   See [synchronize_Q.cxx](../step/initialize_schur_complement_solver/synchronize_Q.cxx).

## Possible improvements

### Better load balancing for BLAS jobs

Currently, each node has `num_primes` BLAS jobs, i.e. `cblas_dsyrk` calls calculating Q_group mod each prime.
These jobs are distributed in a simple round-robin manner: core_1 gets prime_1, core2 - prime_2, etc.

This scheme is certainly not optimal, e.g. if number of cores greatly exceeds number of primes.

For better load balancing, we can split P into M vertical bands P_1..P_M.
Then Q is split into MxM rectangular blocks Q_ij, and one can calculate each `Q_ij mod prime_k` separately.
Since Q is symmetric, it is sufficient to calculate only upper half, `i <= j`.
As a result, we have `num_primes * M * (M + 1) / 2` jobs that we can balance among all the cores on the node.

There are two kinds of jobs:

- Diagonal (square) blocks: `Q_ii = P_i^T P_i` (call `cblas_dsyrk`)
- Off-diagonal blocks: `Q_ij = P_i^T P_j, i < j` (call `cblas_dgemm`)

Generally one can expect that syrk jobs are ~2x faster than gemm, since syrk calculates only the (upper) half of the
matrix.
Then we can assign cost to each job (i.e. n*m for gemm and n(n+1)/2 for syrk, where n and m are block height and width,
respectively) and use any [scheduling algorithm](https://en.wikipedia.org/wiki/Identical-machines_scheduling) to
distribute them. Note that block sizes may differ if width of Q does not divide by M.

Another question is how to choose M.
For example, if we have 100 cores and 99 primes, we are happy with one idle core and do not want to split Q.
If we have 100 cores and 10 primes, then it seems reasonable to set at least M=5 and get 100 BLAS jobs.

Generally, there is a tradeoff between two goals:

1. The jobs should be distributed evenly among the cores.
2. Each BLAS job should be as big as possible (for a single process, one big syrk call is faster that several
   small gemm+syrk calls calculating the same matrix).

If we aim for the second goal only, then M should be set as

```
if num_primes > procsPerNode:
   M = 1
else:
   M = max m: num_primes * m * (m + 1) / 2 <= procsPerNode
M = min(M, Q.Width())
```

If we want better load balancing (at the cost of higher CPU time), then we should use larger M.

Which M is truly optimal, depends on the (unknown) details of `cblas_dsyrk` and `cblas_dgemm` performance for different
matrix sizes.

### RAM optimization: split Q into blocks

For extremely large problems, Q matrix may not fit into a single node.

In this case it also helps to split P into M vertical bands P_1..P_M, so that Q is split into MxM blocks Q_ij.
Then the algorithm can be schematically written as

```
for i,j=1..M:
   - Calculate submatrix Q_group_ij = P_i^T P_j (from all blocks on the node), as in the current algorithm.
   - Send the result to the global Q matrix.
```

With this procedure, memory window size for blocks is reduced by a factor of M.
Memory window size for residues of Q and for Q_group are reduced by a factor of MxM.

Increasing M and the number of cores/nodes, in principle we can fit arbitrarily large problems into memory.