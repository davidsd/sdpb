# Experimental branches and pull requests

*Last updated: 28 November 2026.*

Options and a GPU backend for SDPB, `pmp2sdp` and `spectrum`, offered as pull requests against `master` before they
are reviewed and merged; they may change or be withdrawn. Options are off by default, and with them off the output
is byte-identical to `master`'s (SDPB's nmax = 6 test problem, 1 and 4 ranks) unless the table says otherwise. Every
branch passes SDPB's test suite and adds tests of its own. Until a branch has a human review, compare it with
`master` on a problem whose answer you know before relying on it.

## Credit

You may use these branches **on the condition that you credit the people in the Contributors column of each
branch you use and of every branch it is built on** (its Base, and the Base's bases: each row lists only the people
whose work went into that branch itself), in the acknowledgements or wherever your field
credits software, and cite SDPB as its [Attribution](../Readme.md#attribution) section asks, naming the branch
(e.g. "the `cpu-mm-trsm` branch of SDPB, davidsd/sdpb#NNN").

## Branches

A branch contains its Base, and the Base's bases. Options are explained in `--help` on each branch and in its
pull request, which also says how it was checked.

### `sdpb`: speed and memory on CPUs

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `cpu-mm-trsm` | [#311](https://github.com/davidsd/sdpb/pull/311) | each Schur complement block's Cholesky factor L and the solve L⁻¹B (B the free variable matrix) by exact multimodular arithmetic (integers modulo many primes, multiplied with BLAS). Not bitwise: results differ from `master`'s within a proven error bound | `--multimodularTrsm` | `master` | Vasiliy Dommes (the multimodular method of `bigint_syrk`, which it extends); David Simmons-Duffin | none yet |
| `cpu-mm-cholesky-q` | [#312](https://github.com/davidsd/sdpb/pull/312) | the Cholesky decomposition of Q by multimodular arithmetic. Not bitwise; its accuracy not yet checked independently | `--multimodularCholeskyQ <panel>` | `cpu-mm-trsm` | Vasiliy Dommes (the multimodular method of `bigint_syrk`, which it reuses); David Simmons-Duffin | none yet |
| `cpu-q-solve-options` | [#313](https://github.com/davidsd/sdpb/pull/313) | Q's solves and Cholesky on many nodes: solves on fewer ranks, node-shared diagonal blocks, a solve cache. Bitwise except `--qSolveRanks`, which changes the iterates at working precision | `--qSolveRanks`, `--qCholeskyNodeDiag`, `--qCholeskyNodeFactor`, `--qSolveCache` | `cpu-mm-cholesky-q` | Vasiliy Dommes (the multimodular method of `bigint_syrk`, which it builds on); David Simmons-Duffin | none yet |
| `cpu-q-syrk-paths` | [#314](https://github.com/davidsd/sdpb/pull/314) | faster paths of Q's BigInt syrk on many nodes | `--batchCrt`, `--batchResidues`, `--residuesNoClear`, `--residuesStarVR` | `cpu-q-solve-options` | Vasiliy Dommes (`bigint_syrk`, which it speeds up); David Simmons-Duffin | none yet |
| `cpu-grouped-syrk` | [#329](https://github.com/davidsd/sdpb/pull/329) | Q = PᵀP (P = L⁻¹B) taken by groups of rows of P that have the same nonzero columns, over those columns only: P's exact zeros are skipped, and Q is the same bit for bit. 1.26× per iteration at nmax = 18 on 4 nodes and 1.25× at nmax = 26 on 16, against the same build without the option; none at nmax = 14, where most of P is nonzero. No memory saved | `--groupedSyrk` | `cpu-q-syrk-paths` | Vasiliy Dommes (`bigint_syrk`, which it speeds up); David Simmons-Duffin | none yet |
| `cpu-sparse-b` | [#315](https://github.com/davidsd/sdpb/pull/315) | the free variable matrix B stored by its nonzeros (memory) | `--sparseB` | `cpu-mm-trsm` | David Simmons-Duffin | none yet |
| `lab/pub3-step-length-double` | [#332](https://github.com/davidsd/sdpb/pull/332) | each block's smallest eigenvalue for the step lengths in double, with SDPB's multiprecision eigenvalue only for a block whose double value could move the step length by more than a tolerance. About 7% less time per iteration at 832 bits and 9% at 320 on sdpb-test-04 (2 nodes, with the CPU stack's options on in both arms); 14% on sdpb-test-02 at 320 bits on `master` alone. Not bitwise: the iterates agree with the default's to about 14 digits at 832 bits and 11 at 320 (the accuracy of a double step length), with the same verdict and iteration count in a whole solve. Merge [#331](https://github.com/davidsd/sdpb/pull/331) first: low-precision runs need its z.txt fix | `--stepLengthDouble` (`--stepLengthDoubleTol`, default 1e-10; `--stepLengthCheck` for testing) | `master` | Antoine Tilloy ([SDPviBe](https://github.com/AntoineTilloy/SDPviBe)'s checked float64 step lengths, which the GPU backend took and this branch brings to CPUs: the error bound and the tolerance that decide when double is enough); David Simmons-Duffin | none yet |

### `sdpb`: GPU backend

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `gpu-core` | [#319](https://github.com/davidsd/sdpb/pull/319) | the device library `src/gpu` (CUDA): exact multiprecision products on INT8 tensor cores, block Cholesky | its own build (SDPB's build is untouched) | `master` | Antoine Tilloy ([SDPviBe](https://github.com/AntoineTilloy/SDPviBe), the prototype of the backend, whose algorithms and code it adapts; `docs/GPU_backend_credits.md` on the branch lists what came from where); David Simmons-Duffin (asked for a GPU backend inside SDPB with SDPviBe as its prototype) | none yet |
| `gpu-solver` | [#320](https://github.com/davidsd/sdpb/pull/320) | whole SDPB iterations on one GPU (about 20× per iteration at nmax = 14 on one H200, against SDPB on 5 nodes × 124 ranks), at `--precision` 831 to 862 bits only (`gpu-general-l` adds others) | build with the GPU backend, then `sdpb --gpu` | `gpu-core` | Antoine Tilloy (SDPviBe's solver, which it follows: blocks batched by size, the "lean" mode without a stored P); Zechuan Zheng (the partial-sum matrix-vector products, the sparse storage of B, and the extension of the lean mode to SDPs that fit on the GPU, from his GPU work on SDPviBe); David Simmons-Duffin | none yet |
| `gpu-distributed-q` | [#321](https://github.com/davidsd/sdpb/pull/321) | Q distributed over the GPUs of a node (memory) | `mpirun -n G sdpb --gpu` for G GPUs, with NCCL | `gpu-solver` | Zechuan Zheng (his tiled Q); David Simmons-Duffin | none yet |
| `gpu-general-l` | [#322](https://github.com/davidsd/sdpb/pull/322) | GPU kernels for `--precision` from 704 to 1216 bits | `--gpu --precision ...` | `gpu-distributed-q` | David Simmons-Duffin | none yet |
| `gpu-resident` | [#323](https://github.com/davidsd/sdpb/pull/323) | resident mode on the GPU | `--gpu --controlDir` | `gpu-general-l` and `cpu-resident` | David Simmons-Duffin | none yet |
| `gpu-multinode` | [#324](https://github.com/davidsd/sdpb/pull/324) | Q distributed across nodes | `sdpb --gpu` on several nodes | `gpu-resident` | David Simmons-Duffin | none yet |
| `gpu-fast-q` | [#326](https://github.com/davidsd/sdpb/pull/326) | faster Q on several GPUs, the same results bit for bit: P's residues formed once per rank instead of once per tile, matrix-vector products split over more threads (exact sums), and Q = PᵀP taken by groups of rows with the same nonzero columns. About 1.8× per iteration at nmax = 22 on two H200s, 1.5× on eight; 1.6× at nmax = 18 on two; ~2% on one GPU. Uses ~15 GB more device memory per GPU at nmax = 22 | on by default with `--gpu`; `SDPB_GPU_STRUCTURED_Q=0 SDPB_GPU_Q_SINGLE_PASS=0 SDPB_GPU_GEMV_MIN_TERMS=256 SDPB_GPU_Q_BASE=256` restore the earlier path (see `src/gpu/README.md`) | `gpu-multinode` | Zechuan Zheng (the partial-sum matrix-vector products, which it splits further); David Simmons-Duffin | none yet |
| `gpu-sdpa` | [#336](https://github.com/davidsd/sdpb/pull/336) | SDPA-format input (`.dat-s`: dense, sparse or gzip-compressed, with the existing SDPA grammar and objective-constant convention) solved on the GPU backend, one MPI rank per GPU on a single node. The constraint matrices stay sparse on the device; the Schur complement is assembled from the transformed constraints and handed to the GPU stack's multiprecision factorization, solves and step lengths. Checked on one H100: a 3,410-constraint SDPA problem at 768 bits reaches primal-dual optimality in 73 iterations as an earlier single-GPU SDPA solver does, both objectives agreeing to 25 digits and the solution files byte for byte; about 4.3 s per iteration. Sparse input is parsed on `OMP_NUM_THREADS` threads (dense input still serially). Only one GPU has been tested | `sdpb --gpu -s problem.dat-s` (CPU execution stays the default) | `gpu-fast-q` | Vasiliy Dommes (the SDPA solver and its multimodular Schur construction, which the driver follows); Antoine Tilloy ([SDPviBe](https://github.com/AntoineTilloy/SDPviBe)'s running column maxima and headroom, used for the streamed Gram alignment); Zechuan Zheng (the partial-sum matrix-vector product, adapted for the sparse constraint products, and his streaming Schur construction, which informs releasing the transformed constraints after accumulation); David Simmons-Duffin | none yet |

### `sdpb`: new options

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `sdpb-decision-stop` | [#306](https://github.com/davidsd/sdpb/pull/306) | stop as soon as the iterate decides a bound | `--dualObjectiveStop`, `--primalObjectiveStop` | `master` | David Simmons-Duffin | none yet |
| `initial-seed` | [#304](https://github.com/davidsd/sdpb/pull/304) | a cold start built from the SDP's own rows (fewer iterations on cold calls) | `--initialSeed=rows` | `master` | David Simmons-Duffin; Zechuan Zheng (the finding that equilibrated sample scalings cut SDPB's iteration counts, which this cold start grew out of) | none yet |
| `sdpb-reconstruct-xymu` | [#307](https://github.com/davidsd/sdpb/pull/307) | restart from a text checkpoint that has only x and y | `--reconstructXYmu mu` | `master` | David Simmons-Duffin | none yet |
| `cpu-resident` | [#316](https://github.com/davidsd/sdpb/pull/316) | resident mode: sdpb stays running between the calls of a scan and reads commands from a directory (protocol: `src/sdpb/resident/Readme.md`) | `--controlDir` | `cpu-sparse-b` and `pmp2sdp-reuse-sdp` | David Simmons-Duffin (proposed it: a scan mode inside sdpb, so that it never exits between calls) | none yet |

### `pmp2sdp` and `spectrum`

| Branch | PR | What it does | Turned on by | Base | Contributors | Human review |
|---|---|---|---|---|---|---|
| `pmp2sdp-pmp-info-gatherv` | [#305](https://github.com/davidsd/sdpb/pull/305) | `pmp_info.json` collected exactly, the same on any number of ranks | always on (**changes `pmp_info.json`** written on more than one rank: before, some prefactors were rounded) | `master` | David Simmons-Duffin | none yet |
| `pmp2sdp-reuse-sdp` | [#310](https://github.com/davidsd/sdpb/pull/310) | reuse the blocks of unchanged input files from a previous SDP | `--reuseSdp <sdp>` | `pmp2sdp-pmp-info-gatherv` | David Simmons-Duffin | none yet |
| `pmp2sdp-nodes-at-infinity` | [#317](https://github.com/davidsd/sdpb/pull/317) | interpolation conditions at x = ∞ for reduced prefactors. **Offered for discussion: not shown to be useful** | `--reducedPrefactorNodesAtInfinity M` | `master` | David Simmons-Duffin | none yet |
| `spectrum-zero-test` | [#318](https://github.com/davidsd/sdpb/pull/318) | a second, local test of which minima are zeros | `--zeroTest curvature` (the default stays `geometricMean`) | `master` | Vasiliy Dommes (raised the problem with spectrum's zero test, gave the test cases and suggested approaches); David Simmons-Duffin (chose to keep the default) | none yet |
| `lab/pub4-two-set-encoding` | [#333](https://github.com/davidsd/sdpb/pull/333) | two sets of points: the polynomials stay defined at the sample points, and each block's SDP is written at a second set of encoding points. The polynomial program and its optimum are the same; the SDP is better conditioned, so it can need fewer iterations or less precision. At nmax = 14 (one excluded point of the 3d Ising bootstrap with the stress tensor) 60 iterations at 320 bits against the default's 106 at 832, about 3.2× a whole call (1.7× against the default at 320 bits, where it passes marginally); at nmax = 22 no precision gain (512 bits, as the default), 128 iterations against 184. With the option on, the objective equals the default's to working precision, not byte for byte. With the option off, one change remains, a fix: `spectrum`'s zero test for constant (degree-0) blocks now uses the reduced prefactor at the block's point instead of its stored sample scaling, which removes a spurious zero when a PMP supplies its own scalings (as XML PMPs and SDPB's dfibo test do); PMPs whose constant blocks have scaling 1 are unaffected. Merge [#331](https://github.com/davidsd/sdpb/pull/331) first: at low precision `master` can write a z.txt with n·z far from 1 | `--encodingPoints functionalMatched` (or `chebyshev`, `file`; see `docs/Two_Set_Encoding.md` on the branch) | `master` | David Simmons-Duffin (the idea: one set of points for approximating the blocks, another for encoding the SDP) | none yet |

## Building

Check out a branch once its PR is open ("—" in the PR column: not yet):

    git clone -b <branch> https://github.com/davidsd/sdpb.git

A branch contains its bases, so one checkout gets a whole stack: `gpu-multinode` has every GPU branch (and
`cpu-resident` with its bases), `cpu-grouped-syrk` the multimodular CPU solver. To combine branches from different
stacks, merge them yourself (`git fetch https://github.com/davidsd/sdpb.git <branch> && git merge FETCH_HEAD`).
Many of them add to the same lists (option registration, tests, build files), so expect conflicts there, and
run the tests on the result.

Build with CMake as [Install.md](../Install.md) describes, then run the tests
(`ctest --test-dir build --output-on-failure`).

**The GPU backend** needs NVIDIA GPUs with INT8 tensor cores (developed on H200s, CUDA 12):

    cmake -B build -DSDPB_ENABLE_GPU=ON -DCMAKE_CUDA_ARCHITECTURES=90-real <the usual options>

(`90` is the H100/H200 compute capability; use your GPU's), then run `sdpb --gpu`. Without the option the build
needs no CUDA. `src/gpu/README.md` has the details.

**Runs at low precision** (well below SDPB's usual 768 to 1024 bits) can also need a fix in the bootstrap collaboration's
Elemental fork, which SDPB builds against: its `Axpy` with α = ±1 and its `MakeSymmetric` round away the guard limb, so a
symmetrized matrix differs from its transpose in the last limb, and a Cholesky decomposition can fail (in tests, at 192 bits on
a small problem and at 384 on sdpb-test-01). Merge request
[!11](https://gitlab.com/bootstrapcollaboration/elemental/-/merge_requests/11) of the fork (open)
makes both exact; build Elemental from it for such runs.

## Adding a branch

Other groups' experimental branches are welcome: a pull request that adds one fills in its row of the table.

## Problems

Report problems as an [issue on davidsd/sdpb](https://github.com/davidsd/sdpb/issues) naming the branch, or on its
pull request, with what you ran (command, ranks and nodes, precision) and what you expected.
