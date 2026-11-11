# Resident mode: `sdpb --controlDir <dir>`

One SDPB process loads an SDP once and then executes commands until `quit`: solve from the current iterate, replace
blocks or the objective in memory, set the coefficients of a parametric block, save and restore iterates, write
checkpoints. It is meant for scans that solve many nearby SDPs, e.g. the OPE search (one λ-dependent block, solved
dozens of times per point), where each `sdpb` run would otherwise pay for MPI start-up, reading the SDP and the
checkpoint, writing a checkpoint, and exit. Without `--controlDir`, SDPB solves once as before. Resident mode is for
polynomial matrix programs (pmp2sdp output); with SDPA input `--controlDir` is refused.

## Starting it

```
srun sdpb --sdpDir=<sdp> --checkpointDir=<ck> --precision=960 [any sdpb options] --controlDir=<dir> [--controlIdleTimeout=<s>]
```

Start-up is the usual one: the SDP is read, the initial checkpoint is loaded if there is one (else the initial
point), and the timing run is done if needed. Then SDPB writes `<dir>/ready.json` (host, pid, number of blocks, dual
dimension, precision, ranks) and waits.

## Protocol

- The driver writes the commands `<dir>/cmd.1.json`, `cmd.2.json`, ... Each must appear at once: write it under another
  name in `<dir>`, then rename it.
- SDPB answers each with `<dir>/res.<seq>.json` (also renamed into place): `{"seq": n, "results": [...]}`, one object
  per command, or `{"error": "..."}` for a command that failed.
- Errors: each command first validates its input on all ranks (fields, files, block sizes, N, checkpoint files),
  before it changes anything; such an error is answered and SDPB carries on with its state unchanged. An error later,
  inside the work itself (e.g. a malformed PMP file, a write that fails), may happen on some ranks only, inside a
  collective operation, so SDPB aborts MPI as `sdpb` does; restart from the last checkpoint.
- A command file holds one command (a JSON object with `"cmd"`), or a batch `{"commands": [...]}`, executed in order
  until the first error.
- While waiting, rank 0 polls for the next file (100 µs, doubling to 20 ms) and touches `<dir>/alive` every 10 s. The
  other ranks wait in `MPI_Ibarrier` and sleep, so idle ranks leave the CPUs to the driver.
- At start, SDPB removes `res.*`, `ready.json`, `exit.json`, `stop` and `alive` left in `<dir>` by an earlier
  process; command files are kept (they may be written in advance, as SDPB's integration test does). A driver must not
  leave old `cmd.*` files behind.
- A rank whose parent process (the launcher: slurmstepd, mpirun's daemon) disappears while it waits for a command
  exits, so that a killed launcher leaves no process behind.
- SDPB stops on `quit`, on SIGTERM, when `<dir>/stop` appears, or after `--controlIdleTimeout` seconds without a
  command. In the last three cases it first writes a checkpoint to `--checkpointDir` if the iterate is newer than the
  last one written there (as `sdpb` does on SIGTERM). It writes `<dir>/exit.json` (reason, last seq).

There is no client in this repository; the commands below are the whole protocol.

## Commands

| `cmd` | fields | does |
|---|---|---|
| `solve` | `args`, `inherit`, `start`, `checkpoint`, `write_solution`, `derived` | run the solver |
| `update_pmp` | `pmp` | replace the changed blocks (and the objective) from PMP files, converted in memory |
| `update_sdp` | `sdp`, `blocks` | re-read blocks and the objective from an SDP directory written by pmp2sdp |
| `load_family` | `name`, `block`, `components`, `form_dim`, `normalization_file` | make a block parametric: c = Σ w_i c_i, B = Σ w_i B_i |
| `set_coefficients` | `family`, `coefficients` or `lambda` | set the w_i of a family and rebuild its block |
| `negate_normalization` | | n → −n (e.g. upper vs lower bound of an OPE coefficient at a fixed point) |
| `save_iterate`, `restore_iterate`, `drop_iterate` | `name` or `dir` | keep (x, X, y, Y) in memory, or as a checkpoint in `dir` |
| `checkpoint` | `dir` | write a checkpoint (default: `--checkpointDir`), with `block_timings` |
| `status` | | solves done, iterates, families |
| `quit` | `checkpoint` | stop (optionally after a checkpoint) |

### solve

- `args`: an sdpb command line for this solve, e.g. `["--outDir=out_17", "--maxIterations=500"]`. Options it does not
  give are taken from the resident process's command line. They are parsed at the GMP precision of sdpb's own
  command-line parsing (before `--precision` is set), so BigFloat options (thresholds, centering parameters) have the
  same values as in an `sdpb` run. `--precision`, `--sparseB`, `--maxMemory`, `--maxSharedMemory` and
  `--procGranularity` are fixed at start-up (the BigInt syrk context and its shared memory windows are kept between
  solves). Boolean switches given at start-up cannot be turned off per solve. An `--sdpDir` in `args` must name
  the SDP in memory (the last `update_sdp` directory, or the start-up one); to solve another SDP, update first.
- `inherit` (default true): with `false`, the solve's options are its `args` only, as for an `sdpb` run with that
  command line; nothing is taken from the start-up command line (the options fixed at start-up must still agree, so
  `args` must repeat them). A drop-in replacement for separate `sdpb` runs should use it: with inheritance, a switch
  that the start-up command gave (say `--detectPrimalFeasibleJump`) stays on in every solve. `ready.json` lists
  `"solve_inherit"` in `features`.
- `start`: `"current"` (default: the iterate in memory), `"initial"` (X = Ω_p I, Y = Ω_d I, x = y = 0), `"checkpoint"`
  (load from the solve's `--initialCheckpointDir`/`--checkpointDir`), or `{"iterate": "<name>"}`.
- `checkpoint` (default false): write a checkpoint to `--checkpointDir` after the solve. Without it, no checkpoint is
  written (except on SIGTERM); the iterate stays in memory.
- `write_solution` (default true): write `out.txt`, `y.txt`, `z.txt`, ... to `--outDir`, as `sdpb` does, so a
  driver can read results the way it reads them today. (`iterations.json` and `c_minus_By/` are written in any case,
  as by `sdpb`.)
- `derived`: extra results in the answer. `"y"`, `"z"` (strings at full precision); `"family_values"`: for each family
  and component i, c_i − B_i·y over the block's rows. For a 1x1 constant block this is the functional on the component
  (times the positive sample scaling; `pmp_values` divides it out): for the OPE search, the cut form α·𝒢_ab.
  `write_derived()` in `Resident.cxx` is where more quantities can be added.
- The answer has `terminateReason`, the objectives and errors (strings), `steps`, `outDir`, and timings (`setup`,
  `run`, `checkpoint`, `solution`).
- `maxRuntime` counts from the start of the solve (in `sdpb`, from the start of the process).

### update_pmp (exact)

`pmp`: the new PMP, as for pmp2sdp (a `.json`/`.nsv` file list or a list of files). Like `pmp2sdp --reuseSdp`: files
listed in the current `pmp_info.json` keep their blocks; the others are parsed and converted by pmp2sdp's own code
(`read_polynomial_matrix_program`, `Output_SDP`), on whichever ranks parse them, and each new block goes through the
SDP's block file format (`write_block_data` → `SDP_Block_Data`) to the root of the group that owns it. The objective
is re-read and, if it changed, set through the text of `objectives.json`, as SDPB reads it. So the SDP in memory equals
the one SDPB would read from `pmp2sdp` (`--reuseSdp`) output. Requirements, checked: the same number of blocks, the
reused files at the same block indices, each new block with the old one's size, the same normalization (every block
depends on it). The SDP directory must have `pmp_info.json` (pmp2sdp writes it).

### update_sdp (exact)

`sdp`: a directory written by pmp2sdp (e.g. with `--reuseSdp`). `blocks`: `"changed"` (default: blocks whose source
file in `pmp_info.json` differs from the current one), `"all"`, or a list of indices. Re-reads those blocks and
`objectives.json` as SDPB's constructor does. A changed normalization needs `"all"`.

### load_family, set_coefficients (parametric blocks)

A family makes block `block` parametric: its constraint data become c = Σ_i w_i c_i and B = Σ_i w_i B_i, with the
components (c_i, B_i) converted once from PMP files `components` (one matrix each, with the block's size, sample
points and bilinear basis). The bilinear bases must agree. The components are converted with the SDP's normalization:
give `normalization_file` (the PMP file that holds it) to use the values pmp2sdp parsed, exactly; otherwise the
values of `normalization.json` are used (a decimal round trip). `negate_normalization` negates the families' c_i too;
after `update_sdp` changes the normalization, `set_coefficients` refuses a family loaded before (load it again). With `form_dim` =
m, the components are the upper triangle (row-major) of a symmetric m×m matrix 𝒢, and `set_coefficients` may give
`lambda` (m numbers, decimal or p/q): w_ab = λ_a λ_b, doubled for a < b, so that the block is λᵀ𝒢λ. Otherwise
`coefficients` gives the w_i. The sum is computed in that order at the working precision on the root of the block's
group, then distributed. This is not bitwise equal to contracting λ outside SDPB and converting the result (the
rounding differs: ~1e-194 relative in y at 768 bits in our runs, same verdicts and iteration counts); it is
exact when the w_i select one component.

Clients beyond the OPE scan: any parameter p on which some blocks depend polynomially or linearly, B(p) = Σ f_i(p) B_i
(e.g. the driver computes f_i(p) = p^i), one family per such block.

### negate_normalization (exact)

n → −n. With the same eliminated index, this negates c in every block and the objective constant, and leaves B and b
unchanged; El::BigFloat (GMP mpf) arithmetic is sign-symmetric, so the result equals the SDP pmp2sdp writes for the
negated normalization, bitwise.

### Iterates and checkpoints

`save_iterate {"name"}` copies (x, X, y, Y) in memory (about two X-sized blocks per copy); with `"dir"` it writes a
checkpoint there instead (the main checkpoint's generation counters are kept, as for `checkpoint {"dir"}` with another
directory). `restore_iterate` is the inverse; it first checks that the checkpoint has a file for every rank.
`checkpoint` writes to `--checkpointDir` and copies `block_timings` next to it, so that a new process restarted from
it has the same block distribution.

## Exactness

For the same sequence of SDPs and options, the resident process takes the same steps as `sdpb` runs chained through
checkpoints: the iterate in memory is what the checkpoint would hold (binary checkpoints are exact), the SDP in memory
is what SDPB would read (`update_pmp`, `update_sdp`, `negate_normalization`), the per-solve options have the same
values, and the solver's other state is rebuilt by `SDP_Solver::run()` at each solve, as in a new process. The only
state kept between runs is the BigInt syrk context (shared memory windows), which depends only on the block
structure. Checked outside the test suite: out.txt, y, z and every checkpoint file bitwise equal, call by
call, for a 16-call λ sequence with alternating verdicts (update_pmp and update_sdp), objective changes, and
normalization flips, at 2 ranks; with multi-rank block groups (8 and 24 ranks) and 2 simulated nodes; and for 4
calls of an nmax10 OPE search (113 blocks, 960 bits) on one 128-rank node, where one-step hot
calls took 72.9 s instead of 110.4 s (pmp2sdp --reuseSdp + sdpb).

Replacing a block in place needs `copy_matrix_from_root` to really copy: it used `QueueUpdate`, which adds, so on a
block owned by several ranks the old and new data were summed. A 2-node run showed it (different verdicts from the
second solve on); it now zeros the destination first, which changes nothing when an SDP is read. With the fix the
2-node replay is bitwise for all 16 calls.

The integration test `resident` (`test/src/integration_tests/cases/resident.test.cxx`) solves A, B, A both ways
(SingletScalarAllowed nmax6: excluded, allowed, excluded; the changed block on a multi-rank group) and compares
out.txt, y, z and every checkpoint file after each solve. It fails without the copy fix. It then solves B three more
times from the iterate saved after the first solve (`restore_iterate` by name, `restore_iterate` from a `save_iterate`
checkpoint, and `"start": {"iterate": ...}`), and each must equal the file route's solve of B from the first checkpoint.

## Failure handling (driver side)

- SDPB dies (node failure, MPI error): the `srun` exits; restart from the last checkpoint written (`checkpoint`, or
  `solve` with `"checkpoint": true` every N solves) and resend the current λ. Restarting from the checkpoint of solve k
  and replaying solves k+1, ... gives the same iterates as before (determinism, with the same `block_timings`).
- Time limit: SIGTERM makes SDPB stop the current solve, write a checkpoint and exit; or quit early with
  `{"cmd": "quit", "checkpoint": true}`.
- A command that fails (bad arguments, a file that does not exist) gets an `error` answer; SDPB keeps its state.
  Errors inside collective operations (e.g. one rank failing to read a file) still abort MPI, as in `sdpb`.
