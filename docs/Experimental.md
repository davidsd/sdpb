# Experimental branches and pull requests

*Last updated: 11 November 2026.*

Options and a GPU backend for SDPB, `pmp2sdp` and `spectrum`, offered as pull requests against `master` before they
are reviewed and merged; they may change or be withdrawn. Options are off by default, and with them off the output
is byte-identical to `master`'s (SDPB's nmax = 6 test problem, 1 and 4 ranks) unless the table says otherwise. Every
branch passes SDPB's test suite and adds tests of its own. Until a branch has a human review, compare it with
`master` on a problem whose answer you know before relying on it.

## Credit

You may use these branches **on the condition that you credit the people in the Contributors column of each
branch you use** (the column includes its bases' contributors), in the acknowledgements or wherever your field
credits software, and cite SDPB as its [Attribution](../Readme.md#attribution) section asks, naming the branch
(e.g. "the `cpu-mm-trsm` branch of SDPB, davidsd/sdpb#NNN").

## Branches

A branch contains its Base, and the Base's bases. Options are explained in `--help` on each branch and in its
pull request, which also says how it was checked.

### `sdpb`: speed and memory on CPUs

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `timers-bounded` | [#308](https://github.com/davidsd/sdpb/pull/308) | below `--verbosity debug`, sdpb stops keeping finished timers that nothing reads (memory grows every iteration otherwise) | always on (no output changes) | `master` | David Simmons-Duffin | none yet |
| `cpu-fast-reduce` | [#303](https://github.com/davidsd/sdpb/pull/303) | `restore_and_reduce` visits only each rank's own entries of Q (faster on many nodes) | always on (no output changes) | `master` | Vasiliy Dommes (`restore_and_reduce`, which it changes); David Simmons-Duffin | none yet |
| `cpu-trsm-large` | [#309](https://github.com/davidsd/sdpb/pull/309) | Elemental's `TRSM_LARGE` for the per-block triangular solves (faster when a block's group has a prime number of ranks) | always on (no output changes) | `master` | David Simmons-Duffin | none yet |
| `cpu-mm-trsm` | [#311](https://github.com/davidsd/sdpb/pull/311) | each Schur complement block's Cholesky factor L and the solve L⁻¹B (B the free variable matrix) by exact multimodular arithmetic (integers modulo many primes, multiplied with BLAS). Not bitwise: results differ from `master`'s within a proven error bound | `--multimodularTrsm` | `master` | Vasiliy Dommes (the multimodular method of `bigint_syrk`, which it extends); David Simmons-Duffin | none yet |
| `cpu-mm-cholesky-q` | — | the Cholesky decomposition of Q by multimodular arithmetic. Not bitwise; its accuracy not yet checked independently | `--multimodularCholeskyQ <panel>` | `cpu-mm-trsm` | Vasiliy Dommes (the multimodular method of `bigint_syrk`, which it reuses); David Simmons-Duffin | none yet |
| `cpu-q-solve-options` | — | Q's solves and Cholesky on many nodes: solves on fewer ranks, node-shared diagonal blocks, a solve cache. Bitwise except `--qSolveRanks`, which changes the iterates at working precision | `--qSolveRanks`, `--qCholeskyNodeDiag`, `--qCholeskyNodeFactor`, `--qSolveCache` | `cpu-mm-cholesky-q` | Vasiliy Dommes (the multimodular method of `bigint_syrk`, which it builds on); David Simmons-Duffin | none yet |
| `cpu-q-syrk-paths` | — | faster paths of Q's BigInt syrk on many nodes | `--batchCrt`, `--batchResidues`, `--residuesNoClear`, `--residuesStarVR` | `cpu-q-solve-options` and `cpu-fast-reduce` | Vasiliy Dommes (`bigint_syrk`, which it speeds up); David Simmons-Duffin | none yet |
| `cpu-sparse-b` | — | the free variable matrix B stored by its nonzeros (memory) | `--sparseB` | `cpu-mm-trsm` | Vasiliy Dommes (via `cpu-mm-trsm`); David Simmons-Duffin | none yet |

### `sdpb`: GPU backend

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `gpu-core` | — | the device library `src/gpu` (CUDA): exact multiprecision products on INT8 tensor cores, block Cholesky | its own build (SDPB's build is untouched) | `master` | Antoine Tilloy ([SDPviBe](https://github.com/AntoineTilloy/SDPviBe), the prototype of the backend, whose algorithms and code it adapts; `docs/GPU_backend_credits.md` on the branch lists what came from where); David Simmons-Duffin (asked for a GPU backend inside SDPB with SDPviBe as its prototype) | none yet |
| `gpu-solver` | — | whole SDPB iterations on one GPU (about 20× per iteration at nmax = 14 on one H200, against SDPB on 5 nodes × 124 ranks), at `--precision` 831 to 862 bits only (`gpu-general-l` adds others) | build with the GPU backend, then `sdpb --gpu` | `gpu-core` | Antoine Tilloy (SDPviBe's solver, which it follows: blocks batched by size, the "lean" mode without a stored P; and via `gpu-core`); Zechuan Zheng (the partial-sum matrix-vector products, the sparse storage of B, and the extension of the lean mode to SDPs that fit on the GPU, from his GPU work on SDPviBe); David Simmons-Duffin | none yet |
| `gpu-distributed-q` | — | Q distributed over the GPUs of a node (memory) | `mpirun -n G sdpb --gpu` for G GPUs, with NCCL | `gpu-solver` | Zechuan Zheng (his tiled Q; and via `gpu-solver`); Antoine Tilloy (via `gpu-solver`); David Simmons-Duffin | none yet |
| `gpu-general-l` | — | GPU kernels for `--precision` from 704 to 1216 bits | `--gpu --precision ...` | `gpu-distributed-q` | Antoine Tilloy, Zechuan Zheng (via `gpu-distributed-q`); David Simmons-Duffin | none yet |
| `gpu-resident` | — | resident mode on the GPU | `--gpu --controlDir` | `gpu-general-l` and `cpu-resident` | Antoine Tilloy, Zechuan Zheng (via `gpu-general-l`); Vasiliy Dommes (via `cpu-resident`); David Simmons-Duffin | none yet |
| `gpu-multinode` | — | Q distributed across nodes | `sdpb --gpu` on several nodes | `gpu-resident` | Antoine Tilloy, Zechuan Zheng, Vasiliy Dommes (via `gpu-resident`); David Simmons-Duffin | none yet |

### `sdpb`: new options

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `sdpb-decision-stop` | [#306](https://github.com/davidsd/sdpb/pull/306) | stop as soon as the iterate decides a bound | `--dualObjectiveStop`, `--primalObjectiveStop` | `master` | David Simmons-Duffin | none yet |
| `initial-seed` | [#304](https://github.com/davidsd/sdpb/pull/304) | a cold start built from the SDP's own rows (fewer iterations on cold calls) | `--initialSeed=rows` | `master` | David Simmons-Duffin; Zechuan Zheng (the finding that equilibrated sample scalings cut SDPB's iteration counts, which this cold start grew out of) | none yet |
| `sdpb-reconstruct-xymu` | [#307](https://github.com/davidsd/sdpb/pull/307) | restart from a text checkpoint that has only x and y | `--reconstructXYmu mu` | `master` | David Simmons-Duffin | none yet |
| `cpu-resident` | — | resident mode: sdpb stays running between the calls of a scan and reads commands from a directory (protocol: `src/sdpb/resident/Readme.md`) | `--controlDir` | `cpu-sparse-b` and `pmp2sdp-reuse-sdp` | David Simmons-Duffin (proposed it: a scan mode inside sdpb, so that it never exits between calls); Vasiliy Dommes (via `cpu-sparse-b`) | none yet |

### `pmp2sdp` and `spectrum`

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `pmp2sdp-pmp-info-gatherv` | [#305](https://github.com/davidsd/sdpb/pull/305) | `pmp_info.json` collected exactly, the same on any number of ranks | always on (**changes `pmp_info.json`** written on more than one rank: before, some prefactors were rounded) | `master` | David Simmons-Duffin | none yet |
| `pmp2sdp-reuse-sdp` | [#310](https://github.com/davidsd/sdpb/pull/310) | reuse the blocks of unchanged input files from a previous SDP | `--reuseSdp <sdp>` | `pmp2sdp-pmp-info-gatherv` | David Simmons-Duffin | none yet |
| `pmp2sdp-nodes-at-infinity` | — | interpolation conditions at x = ∞ for reduced prefactors. **Offered for discussion: not shown to be useful** | `--reducedPrefactorNodesAtInfinity M` | `master` | David Simmons-Duffin | none yet |
| `spectrum-zero-test` | — | a second, local test of which minima are zeros | `--zeroTest curvature` (the default stays `geometricMean`) | `master` | Vasiliy Dommes (raised the problem with spectrum's zero test, gave the test cases and suggested approaches); David Simmons-Duffin (chose to keep the default) | none yet |

## Building

Check out a branch once its PR is open ("—" in the PR column: not yet):

    git clone -b <branch> https://github.com/davidsd/sdpb.git

A branch contains its bases, so one checkout gets a whole stack: `gpu-multinode` has every GPU branch (and
`cpu-resident` with its bases), `cpu-q-syrk-paths` the multimodular CPU solver. To combine branches from different
stacks, merge them yourself (`git fetch https://github.com/davidsd/sdpb.git <branch> && git merge FETCH_HEAD`).
Many of them add to the same lists (option registration, tests, build files), so expect conflicts there, and
run the tests on the result.

Build with CMake as [Install.md](../Install.md) describes, then run the tests
(`ctest --test-dir build --output-on-failure`).

**The GPU backend** needs NVIDIA GPUs with INT8 tensor cores (developed on H200s, CUDA 12):

    cmake -B build -DSDPB_ENABLE_GPU=ON -DCMAKE_CUDA_ARCHITECTURES=90-real <the usual options>

(`90` is the H100/H200 compute capability; use your GPU's), then run `sdpb --gpu`. Without the option the build
needs no CUDA. `src/gpu/README.md` has the details.

## Adding a branch

Other groups' experimental branches are welcome: a pull request that adds one fills in its row of the table.

## Problems

Report problems as an [issue on davidsd/sdpb](https://github.com/davidsd/sdpb/issues) naming the branch, or on its
pull request, with what you ran (command, ranks and nodes, precision) and what you expected.
