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
(see [reduce_scatter_DistMatrix.hxx](reduce_scatter_DistMatrix.hxx), previously called `synchronize_Q`).

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

5. All cores are working in parallel to square each matrix and store the results
   in [Residue_Matrices_Window](Residue_Matrices_Window.hxx):

```
Q_group mod prime_1
Q_group mod prime_2
Q_group mod prime_3
...
```

Each `Q_group mod prime_k` is calculated either via single `cblas_dsyrk()` call or
several `cblas_dsyrk()`/`cblas_dgemm()` calls, if Q is split into blocks for better parallelization (see below).

6. `Q_group` is now a `DistMatrix<BigFloat>` distributed over all cores of a node. If some element `Q_group(i,j)` is
   stored e.g. on core 1, this core restores its value from the `Residue_Matrices_Window` using CRT.
7. Restore Q_group (see `Matrix_Normalizer.restore_Q`), i.e. divide by 2^2N and remove normalization.
8. Calculate global Q, which is `DistMatrix<BigFloat>` distributed over all cores, as a sum of all Q_groups.
   See [reduce_scatter_DistMatrix.hxx](reduce_scatter_DistMatrix.hxx).

### Distributing BLAS jobs

How can we distribute BLAS jobs from step 6 among the processes?

The simplest way is to assign each prime to some core
in a simple round-robin manner: core_1 gets prime_1, core2 - prime_2, etc.
Then each core calls `cblas_dsyrk` calculate Q_group modulo given prime.
This scheme is certainly far from optimal, e.g. if the number of cores greatly exceeds the number of primes.

We adopt the following algorithm:

1. Distribute uniformly across all the cores as many primes as we can.
   For each prime p_k, Q_group mod p_k is calculated via single `cblas_dsyrk` call.
2. For each of the remaining `num_primes % num_ranks` primes,
   we can split P into M vertical bands P_1..P_M.
Then Q is split into MxM rectangular blocks Q_ij, and one can calculate each `Q_ij mod prime_k` separately.
Since Q is symmetric, it is sufficient to calculate only upper half, `i <= j`.
   As a result, we have `M * (M + 1) / 2` jobs for each of the remaining primes.

There are two kinds of jobs:

- Diagonal (square) blocks: `Q_ii = P_i^T P_i` (call `cblas_dsyrk`)
- Off-diagonal blocks: `Q_ij = P_i^T P_j, i < j` (call `cblas_dgemm`)

We estimate the cost of a job as a number of matrix elements that it calculates,
i.e. `n*m` for gemm and `n*(n+1)/2` for syrk, where `n` and `m` are block height and width, respectively. This estimate
is correct for naive matrix multiplication and does not account for specific BLAS optimization.
Note that syrk cost is n(n+1)/2, since it calculates only the upper half of the matrix.

3. How do we choose the split factor M?
   The minimal reasonable value can be determined as

```
   min_split_factor = max m: num_jobs(m) <= num_ranks; m = 1..Q.Height()
```

where `num_jobs(m) = (num_primes % num_ranks) * m * (m+1) / 2` is the number of jobs for the remaining primes.
It ensures that each core gets no more than one extra job, thus avoiding extra overhead.

Then, for each M from `min_split_factor` to `min_split_factor+4`,
we create job schedule using
greedy [Longest-processing-time-first scheduling algorithm](https://en.wikipedia.org/wiki/Longest-processing-time-first_scheduling),
(see [LPT_scheduling.hxx](blas_jobs/LPT_scheduling.hxx))
and choose the best schedule.
The best schedule is the one that has minimal total execution time, i.e.,
minimal [max_cost()](blas_jobs/Blas_Job_Schedule.hxx) among all ranks.

In most cases, testing the lowest five M gives almost uniform (up to several percent) load balancing.
Of course, one could get better balancing by increasing M and splitting Q into smaller blocks (1x1 in the limiting
case).
But it would likely increase the total execution time because one big BLAS call is faster than several smaller ones.

## Possible improvements

### Better job scheduling

If benchmarks show that job balancing is bad for some realistic cases, we can consider more sophisticated load
balancing - e.g. use different splitting algorithm or more realistic BLAS job costs (although they are
implementation-specific).

### Block timings

The old way
of [distributing SDP blocks among the cores](../../../Block_Info/allocate_blocks/compute_block_grid_mapping.cxx)
relies on [block timings](../../../../sdpb/write_timing.cxx), which include calculating contribution to Q from a single
block, `run.step.initializeSchurComplementSolver.Q.syrk_XXX`.
Now all blocks from a node are processed together.
So we may include some timing estimate for this step,
e.g. `(Q.syrk for a node) * (block size) / (total size for all blocks on a node)`.

### RAM optimization: split Q into blocks

For extremely large problems, Q matrix may not fit into a single node.

In this case, it also helps to split P into M vertical bands P_1..P_M, so that Q is split into MxM blocks Q_ij.
Then the algorithm can be schematically written as

```
for i,j=1..M:
   - Calculate submatrix Q_group_ij = P_i^T P_j (from all blocks on the node), as in the current algorithm.
   - Send the result to the global Q matrix.
```

With this procedure, memory window size for blocks is reduced by a factor of M.
Memory window size for residues of Q and for Q_group are reduced by a factor of MxM.

Increasing M and the number of cores/nodes, in principle, we can fit arbitrarily large problems into memory.