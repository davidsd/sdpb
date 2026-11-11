# Options for large SDPs on many nodes (experimental)

These options speed up the stages of an SDPB iteration that work on the global N x N matrix Q of the Schur complement
solver after it is computed: its Cholesky decomposition (`Cholesky_Q`) and the two triangular solves with its factor in
each search direction. On many nodes these stages are latency- and memory-bound and get slower as ranks are added. All
options are off by default; with them off, SDPB runs its default code paths. Each one is checked by unit tests and by
bit-for-bit comparisons of whole SDPB runs (iterations, y, x) on 1-16 Expanse nodes.

Measured on Expanse, 128 ranks per node (`srun`), 832 bits, synthetic nmax22 SDP (N = 7,990) on 16 nodes
unless stated; all with `--multimodularTrsm`. The numbers come from single allocations and move by ~5-10% between
allocations.

A combination for 10-16 nodes:

    sdpb ... --multimodularTrsm --multimodularCholeskyQ 64 --qSolveRanks 64 \
      --qCholeskyNodeDiag --qCholeskyNodeFactor --qSolveCache

## `--multimodularCholeskyQ <panel width>`

The Cholesky decomposition of Q as `El::cholesky::UpperVariant3Blocked` (the same panels and redistributions), with the
trailing update `A22 -= A12^T A12` and the diagonal blocks by multimodular BLAS
(`src/sdp_solve/SDP_Solver/run/multimodular/multimodular_cholesky_dist.cxx`). Panel 64 was fastest on 5-16 nodes.
Changes the factor at working precision (like El's own change of block size). nmax22, 16 nodes: 73.0 -> 58.5 s/iter.

## `--qSolveRanks <s>`

The triangular solves with Q's factor (four per iteration) on the first s ranks, in node order (s = 32 or 64: part
of one node; `src/sdp_solve/SDP_Solver/run/step/Q_Grids.hxx`). The factor is moved there by one AllToAll of its upper
triangle; dy is summed over all ranks by one AllReduce and broadcast back. On the full grid each of El's Trsv blocks
gathers a 128 x 128 diagonal block to every rank, which dominates at many nodes. It changes the iterates at working
precision: the solve grid's shape and the AllReduce of dy set the order of the partial sums. nmax22, 16 nodes, with
`--multimodularCholeskyQ 64`: Q solve 12.1 -> 3.5 s, 58.5 -> 50.0 s/iter (`--qSolveRanks 32`). Memory: the s ranks
hold the whole factor (N^2 BigFloats in all) and receive buffers for half of it, on top of Q on all ranks; this is
counted in SDPB's memory estimate (`--maxMemory`).

## `--qCholeskyNodeDiag`

With `--multimodularCholeskyQ`: each panel's diagonal block is copied to every rank once per node instead of once per
rank: every rank serializes its entries into node-shared memory, one leader per node exchanges the nodes' parts, and
every rank reads the block from there (`multimodular/Node_Diag_Exchange.hxx`). Same bits. Pays only across nodes:
-1.8 s per iteration at nmax22 on 16 nodes; -0.4 to -0.9 s per Cholesky_Q at 10 nodes (N = 5,458); nothing on one
node (its extra node barrier costs a little there).

## `--qCholeskyNodeFactor`

With `--multimodularCholeskyQ`: the factorization of each 64 x 64 diagonal block, which every rank repeats, is split
among the node's ranks (the triangular solve of its top-level split by columns, its multimodular GEMM by rows) and
exchanged through node memory. Same bits. Halves that factorization: -0.15 to -0.3 s per Cholesky_Q at N = 5,458 on
one node, -0.5 s per iteration at nmax22 on 16 nodes.

## `--qSolveCache`

The Q solves are El's `cholesky::SolveAfter(UPPER, NORMAL)` line for line, except that the copy of each 128-row
diagonal block of the factor to every rank of the solve grid is made once per iteration (once per node) instead of
once per triangular solve (`step/Q_Solve_Cache.hxx`). Same bits. Best with `--qSolveRanks`: it moves the best solve
grid from 32 to 64 ranks. Costs N x 128 BigFloats per solve rank (~150 MB at nmax22).

Together, `--qCholeskyNodeDiag --qCholeskyNodeFactor --qSolveCache` with `--qSolveRanks 64` instead of 32: nmax22, 16
nodes, 50.0 -> 45.9 s/iter. With `--qSolveRanks 32` on a real nmax18 SDP (3d Ising with the stress tensor, an OPE scan
point, 10 nodes): 33.9 -> 32.3 s/iter, bit for bit the same iterates.

## Testing on one machine

`SDPB_FAKE_NODES=k` (environment, test hook) splits the ranks of each machine into k contiguous groups that SDPB treats
as nodes, each with its share of the machine's memory, so that the multi-node code paths (the reduce over nodes in
`restore_and_reduce`, the per-node copies above) can be tested on one host. The groups never span machines. An
end-to-end test uses it; the unit tests of the node options split the ranks through
`node_diag_exchange_test_node()` instead.
