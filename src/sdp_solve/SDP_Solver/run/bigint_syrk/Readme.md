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

Class [Fmpz_Comb](fmpz/Fmpz_Comb.hxx), which is a wrapper of FLINT's `fmpz_comb_t`,
helps us with steps 2, 3, and 5.

## Going parallel

In SDPB, everything is more complicated than described above, because P is not a `Matrix<BigFloat>`, but
a `vector<DistMatrix<BigFloat>>`.

This means that P is split into horizontal bands (or blocks), and each block is stored on one or several cores.

In turn, Q is `DistMatrix<BigFloat>` distributed over all cores.

In the old SDPB algorithm, each core (or group of cores) was calculating a contribution to Q (
called `Q_group`) from its blocks. Adding up all Q_groups, one got the resulting Q
(see [synchronize_Q.cxx](https://github.com/davidsd/sdpb/blob/3.0.0/src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/synchronize_Q.cxx)).

This is inefficient in terms of memory, since we allocate memory for `Q_group` for each group of cores.

Now we can avoid this duplication, if all cores on a single node are working together on the blocks,
using [MPI shared memory window](https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPI_Win_allocate_shared.3.html).

The parallel algorithm goes as follows (see [compute_Q](../step/initialize_schur_complement_solver/compute_Q.cxx)
and [bigint_syrk_blas](BigInt_Shared_Memory_Syrk_Context/bigint_syrk_blas.cxx) functions):

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

6. `Q_group` is now stored implicitly, as a residues in the `Residue_Matrices_Window`. If some rank on a node needs some element `Q_group(i,j)`, it can restore it from the residues using CRT.
7. Reduce-scatter: calculate global Q, which is a `DistMatrix<BigFloat>` distributed over all cores, as a sum of all
   Q_groups.
   See [restore_and_reduce.cxx](BigInt_Shared_Memory_Syrk_Context/restore_and_reduce.cxx).
8. Restore Q (see `Matrix_Normalizer.restore_Q`), i.e. divide by 2^2N and remove normalization.

Steps 2, 3 and 8 are performed in [compute_Q()](../step/initialize_schur_complement_solver/compute_Q.cxx) function,
steps 4-7 - in [BigInt_Shared_Memory_Syrk_Context::bigint_syrk_blas()](BigInt_Shared_Memory_Syrk_Context/bigint_syrk_blas.cxx) function.

## Performance and memory optimizations

### Distributing BLAS jobs

How can we distribute BLAS jobs from step 5 among the processes?

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
i.e. `n*m` for gemm and `n*(n+1)/2` for syrk, where `n` and `m` are block height and width, respectively.
This estimate is correct for naive matrix multiplication and does not account for specific BLAS optimizations.
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

### Optimizing shared memory access time

Accessing shared memory window can be significantly slower that accessing local memory (~10x for Expanse HPC).
On systems with NUMA architecture, access time also varies for different cores and different memory regions. It depends
on memory pinning (a memory page gets mapped into the local memory of the processor that first touches it).
For more details, see
e.g. [A Performance Evaluation of MPI Shared Memory Programming, David Karlbom (2016)](https://www.diva-portal.org/smash/get/diva2:938293/FULLTEXT01.pdf).

Two optimizations:

1. Make a "first touch" of memory window, so that each memory region will be pinned to a core that uses it most.
2. Bulk read/write to minimize memory access overhead.

#### First touch

When making a first touch for input window, we can optimize it either for writing block residues or for BLAS calls.

Memory access pattern for computing block residues is quite mosaic-like.
Each process computes residues for all block elements that it owns, and writes them to the shared memory window. Only
residues of a single column modulo single prime are stored consecutively (if the process owns the whole column). In
total, each process is accessing at least to `block_width * num_primes` disconnected memory regions.

For BLAS calls, memory access is more regular:
matrix of residues of P modulo given prime is processed either as a whole, or split into several submatrices.

In our benchmarks, BLAS calls took more time than computing residues.
Thus, in `BigInt_Shared_Memory_Syrk_Context` we make the first touch according to BLAS job schedule. We do the same for
the output memory window.

P.S. Note that memory is pinned page by page, where page size is usually 4096B = 512 doubles, so the memory access is
still non-optimal, especially for smaller blocks.

P.P.S. Creating and touching shared memory window introduces overhead for SDP solver start. For example, initializing
70GB window on Expanse HPC takes about 30 seconds.
This can be significant e.g. for Skydiving algorithm, where solver is restarted after several iterations.

#### Bulk memory access

If we write residues to the shared memory window one by one, it can be rather slow (even slower than BLAS calls).

To optimize it, we compute residues for contiguous memory area (i.e. residues module single prime for each column of a
single block or several consecutive blocks) locally, and then `memcpy` them to the window.
See [compute_block_residues.cxx](BigInt_Shared_Memory_Syrk_Context/compute_block_residues.cxx).

**TODO**: currently, each block are stored in `DistMatrix<BigFloat, MC, MR>`, which means that both columns and rows are
distributed among 2D process grid.
If block is owned at least by four processes, then elements of a column are also distributed among several processes. In
this case, there aren't any contiguous memory areas at all, and the `memcpy` trick is useless.

Possible solutions:

1. Store blocks in `DistMatrix<BigFloat, STAR, VR>`, which means that each column is always owned by a single rank. See
   Fig. 11 and Fig.12 (page 21) of https://www.cs.utexas.edu/~flame/pubs/Elemental1.pdf. NB: this may affect e.g.
   Cholesky decomposition performance
   in [initialize_schur_off_diagonal()](../step/initialize_schur_complement_solver/compute_Q.cxx).
2. Reorder rows of P matrix to ensure that rows coming from each rank are consecutive. NB: it can be non-trivial (and
   cause subtle bugs) when some blocks are owned by a single ranks and others are distributed.

### Reducing memory usage

Storing (double-precision) residues of a BigFloat matrix requires ~4x more memory than the BigFloat matrix itself.
For large problems, shared memory windows for P and/or Q do not fit into memory. Note also that sometimes HPC
administrators may set maximal window size e.g. to 50% of available RAM.

If Q is small and P is big, one can circumvent this issue by increasing number of nodes. However, this introduces extra
communication overhead.
If Q is too big, then one has to split it anyway.

#### Splitting P

If all residues do not fit into memory, then we choose split factor `s` and allow for each MPI group `g` to write only `H_g / s` rows at once to the window, where `H_g` is the total height of all blocks owned by this group.

Then the matrix multiplication [algorithm](#going-parallel) should be modified as:

```
1. Choose a set of primes.
2. Calculate column norms of P.
3. Normalize P and multiply by 2^N
4. Fill Q window for with zeros.
5. While not all rows of P are processed:
   5.1 Each MPI group takes top (H_g/s) remaining rows from its blocks and writes their residues to the P window.
   5.2 Call BLAS jobs to update the Q window.
6. Compute Q_group from the residues stored in the Q window.
7. Reduce-scatter Q_group from all nodes to the global Q.
8. Restore Q (divide by 2^2N and remove normalization).
```

NB: distributed blocks require extra care, to ensure that we're not mixing elements from different rows.

Splitting P saves memory, but introduces some extra overhead:

- More synchronization points: before calling BLAS jobs, we have to wait until all ranks fill input window.
- [Bulk memory access](#bulk-memory-access) is less effective since we're copying small chunks of data each time.
- We make more BLAS calls (for smaller matrices), which should be less effective than fewer BLAS calls (for bigger
  matrices).

P.S. One could also split P vertically and store residues only for several columns.
But this would require to calculate residues for each element multiple times.

#### Splitting Q

Imagine that we split P into M vertical bands P_1..P_M, so that Q is split into MxM blocks Q_ij.
Then the [algorithm](#going-parallel) can be schematically written as

```
for i,j=1..M:
   - Calculate submatrix Q_group_ij = P_i^T P_j (from all blocks on the node), as described above.
   - Reduce-scatter Q_group_ij from all nodes to the global Q_ij.
```

Here we don't need to allocate memory for the whole Q_group and for its residues,
since we can reuse the same buffer for each Q_ij.

If we are splitting both P and Q memory windows, then the algorithm reads:

```
1. Choose a set of primes.
2. Calculate column norms of P.
3. Normalize P and multiply by 2^N.
4. For each i,j=1..M:
   4.1 Fill Q window for with zeros.
   4.2 While not all rows of P are processed:
     4.2.1 Each MPI group takes top (H_g/s) remaining rows from its blocks and writes their residues (for corresponding columns) to the P window.
     4.2.2 Call BLAS jobs to update Q window.
   4.3 Compute Q_group_ij from the residues stored in Q window.
   4.4 Reduce-scatter Q_group_ij from all nodes to the global Q_ij.
5. Restore Q (divide by 2^2N and remove normalization).
```

#### How to choose split factors for P and Q

1. Determine memory limit for shared memory windows. Currently, it is set by command-line option, e.g. `--maxSharedMemory=128G`. If the limit is set to `0` or not provided, SDPB will determine it automatically based on expected memory usage and amount of available RAM. 
2. Choose minimal split factor for Q so that it fits into the memory limit (and leaves some room for P, namely 1 row per each MPI group).
3. Choose minimal split factor for P so that it fits into remaining memory.

Note that our top priority is to minimize split factor for Q,
because splitting Q introduces more significant overhead: (1) expensive reduce-scatter calls and (2) residues for each P
element being calculated several times.

Note also that split factor of Q should be the same for all nodes, since we perform global reduce-scatter for each submatrix of Q.
At the same time, split factor for P differe among the nodes (but should be the same for all ranks on a node).

### Block distribution

SDP blocks are [distributed among the cores](../../../Block_Info/allocate_blocks/compute_block_grid_mapping.cxx)
according to [block timings](../../../../sdpb/write_timing.cxx).
In the old algorithm, for each block XXX, its cost is calculated as a sum
of three timers:

```
sdpb.solve.run.iter_2.step.initializeSchurComplementSolver.Q.cholesky_XXX
sdpb.solve.run.iter_2.step.initializeSchurComplementSolver.Q.solve_XXX
sdpb.solve.run.iter_2.step.initializeSchurComplementSolver.Q.syrk_XXX
```

In the new algorithm, all blocks from a node are processed together, so we cannot measure `Q.syrk_XXX` for a single
block.
Now we use the same two timers `Q.cholesky_XXX` and `Q.solve_XXX`
as before. Instead of `Q.syrk_XXX`, we take time for computing residues and for calling BLAS, and split it among all
blocks proportionally to the block sizes.

**TODO:**
There are other steps in the algorithm which are performed for each block separately.
Ideally we should account for all of them (NB: excluding waiting time!).

**TODO:**
Current block distribution algorithm does not require uniform memory distribution explicitly. Generally, block timings
correlate with block sizes, but two block of the same size can have different timings if one of them contains lots of
zeros.Ideally, we should aim for both optimal memory distribution (to fit the problem into as few nodes as possible) and
timing distribution (to minimize computation time).
