# Real-Q accuracy regression

The original inverse-product Cholesky loses accuracy on the frozen bootstrap
Q matrices. The default corrected mode uses exact whole-matrix power-of-two
scaling and one panel residual correction. This directory adds an independent
MPFR gate and all-entry GPU residual checks, with independent GMP dots.  The frozen real dimensions 3406/5458 pass at both 1024/1280-bit MPFR
reference precisions. The reference columns agree to about 760/675 bits;
all first-32 and 128 spread candidate columns are better than SDPviBe.
This is factor/inverse validation, not a production solver integration.
Full-precision correction is the default. Adaptive correction precision,
exact omission of triangular zeros, and packed symmetric updates are explicit
options. `Factor_Options{false, 0}` selects
the original unrefined recursion for reproducible comparisons.

The optimized native setting is `Factor_Options{true,1,0,-1,1024,true}`
with constructor digits=112 (batch=1 for packed updates). Its worst errors
on 128 spread inverse columns are 2^-629.67 / 2^-545.66, versus SDPviBe
2^-590.50 / 2^-505.27. Direct extraction from stored X gives the same
reported errors as X B with identity columns. Main digits=108 also passes
these real fixtures and all six stress fixtures, with smaller margins.
The corresponding command-line additions are `--correction-digits -1
--skip-triangular-zeros 1 --gram-min 1024 --digits 112` (108 for the reduced
option). The timing harness also uses `--retain-pool 1`.

On one H200, matched native speed ratios are about 1.58x/1.46x SDPviBe;
108 digits gives about 1.66x/1.55x. The original speed ratios were not
retained. Full-correction defaults add 20--23% to the already-accurate real
batched Schur workloads, while still beating SDPviBe about 1.94x/2.24x.
The explicit old option can preserve checked Schur behavior; the real Q
fixtures require correction. No automatic policy is implemented.

The unsuccessful recursive-panel, tiled-update, and low-precision inverse
prototypes are not included. Their command-line flags are rejected instead of
being silently ignored.

SDPviBe by Antoine Tilloy (MIT) is the prototype. The reference baseline is its
unchanged `0d02884` revision. The scalar MPFR Cholesky/forward substitution in
`ref_solve.c` is an independent implementation, extended with checked column
caches and executable comparison gates.

## Real fixtures (not distributed)

The accuracy figures above were measured on two Q matrices taken from SDPB
bootstrap runs, of dimensions 3406 and 5458 (padded internally to 3424 and
5472). These matrices are not distributed with SDPB. The checks that read them
(`q_check`, `real_block_check`, `residual_check` and `ref_solve` on real data)
need a fixture directory in the format below; the synthetic SPD fixtures of the
next section exercise the same checks without them. The checksums of the
fixture files are in `real-q-fixtures.sha256`. Each fixture directory
contains:

- `S0.m.bin`, `S0.e.bin`: symmetric Q, row-major, 28 signed int32 limb planes
  and one int64 exponent per entry.
- `B0.m.bin`, `B0.e.bin`: first 32 identity columns, row-major.
- `P0.m.bin`, `P0.e.bin`: unchanged SDPviBe's L^-1 B, same layout.
- `meta.json`: dimensions and provenance.

A value is the two's-complement integer mantissa times 2^exponent; zero uses
exponent -2^40. Factor and inverse exports from `q_check` are instead padded
**column-major** matrices. The two layouts must not be interchanged. The input
matrices are never shifted or modified.

## Build and run

From `src/gpu`, build the GPU checks with the usual `GMP_DIR` and `ARCH=90`.
`make mpfr-checks MPFR_DIR=<prefix>` builds the two CPU tools separately. Run
large MPFR reference computations on a compute node, not a login node. GPU
commands need a GPU. Set `Q_DIR`, `Q_N`, and an existing `RUN_DIR`.

```
build/q_check --dir "$Q_DIR" --nb 1 --P "$Q_N" --N 32 --K 32 --W 0 \
  --equilibrate 1 --panel-refine 1 --correction-digits -1 --reps 5 \
  --out "$RUN_DIR/cols.bin" --factor-out "$RUN_DIR/L" --inverse-out "$RUN_DIR/X"
build/residual_check --dir "$Q_DIR" --P "$Q_N" --factor "$RUN_DIR/L" --max-log2 -850
build/residual_check --dir "$Q_DIR" --P "$Q_N" --factor "$RUN_DIR/L" \
  --inverse "$RUN_DIR/X" --max-log2 -450
build/ref_solve "$Q_DIR" 0 "$RUN_DIR/cols.bin" 1280 \
  --save-ref "$RUN_DIR/ref1280.txt" --require-no-worse-sv 1
build/ref_solve "$Q_DIR" 0 "$RUN_DIR/cols.bin" 1024 \
  --save-ref "$RUN_DIR/ref1024.txt" --require-no-worse-sv 1
build/compare_references "$Q_N" "$RUN_DIR/ref1024.txt" "$RUN_DIR/ref1280.txt"
```

The MPFR gate checks every selected candidate column against the independently
computed answer and fails if any is worse than SDPviBe. `compare_references`
prints convergence column by column; its optional fourth argument enforces a
log2 error threshold. Choose the reference precision so this convergence is
well beyond the candidate error. A cache is valid only for exactly the same
Q and B; column IDs may be permuted. Missing or duplicate cache columns fail.

The factor negative control must fail (exit 1):

```
build/residual_check --dir "$Q_DIR" --P "$Q_N" --factor "$RUN_DIR/L" \
  --perturb-bits 760 --max-log2 -850
```

`--max-log2` thresholds in these residual commands are regression floors, not
claims of an exact all-entry maximum. The GPU evaluator itself rounds. It is
supplemented by 2048-bit GMP dots at the worst reduction blocks and seeded
positions.
`real_block_check` retains its small sampled smoke test; set its cross-code
agreement threshold to the accuracy of the frozen SDPviBe answer, not to an
unachievable fixed -700 on every condition number. MPFR is the forward gate.

## Exact random stress fixtures

```
python3 test/make_spd_fixtures.py "$RUN_DIR/synthetic"
```

This generates dimensions 64,128,256 with a 200-bit column scale span and a
0- or 100-bit row span. Gram formation uses exact integers, all values fit the
stored precision, and modular nonsingularity proves SPD. No ridge hides the
small eigenvalues. For each fixture, `bench/sdpvibe_q_bench.py --dir <fixture>
--extended-dir <fixture>/reference --columns 32 --sdpvibe <unchanged clone>
--out <fixture>/sv.L` generates B and the independent SDPviBe P. Then use the
same GPU and MPFR commands as for Q. `--columns 128` selects 128 spread identity
columns for the large real matrices, including the last column.

The known-integer and permuted-cache controls run with
`python3 test/check_reference_cache.py build/ref_solve <temporary directory>`.
Each gate above has a negative control: a deliberately perturbed input or
reference that the gate must reject.
