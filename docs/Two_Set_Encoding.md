# Two-set encoding in `pmp2sdp` (`--encodingPoints`)

*Experimental (offered as a pull request, off by default). With the options off, `pmp2sdp` writes the same files as
`master`, byte for byte.*

## What it does

`pmp2sdp` turns each positive matrix with prefactor
$\chi(x)\, \sum_n y_n P_n(x) \succeq 0$ for $x \ge 0$ into an SDP block, by sampling at $d+1$ points, where $d$ is the
block's degree. By default one set of points does two jobs:

1. **It defines the polynomial program.** With a `reducedPrefactor` (given in the input, or made by `--maxNumPoles`),
   the block's polynomials are replaced by the degree-$d$ polynomial $R(x)$ that interpolates
   $P(x) / \prod(\text{removed poles})$ at the sample points (arXiv:2509.14307, Section 4.3). These are the
   **approximation points**: the choice of these points sets how well the reduced program approximates the original
   one.
2. **It encodes the SDP.** The rows of $B$ and $c$ are $\sigma_k R(x_k)$ at the same points $x_k$, with the sample
   scalings $\sigma_k$, and the bilinear bases are polynomials orthonormal for the discrete measure
   $\sum_k \sigma_k \delta(x - x_k)$ (and $\sum_k \sigma_k x_k \delta(x - x_k)$ for the odd part). The choice of these
   points and scalings sets how well conditioned the SDP is, but not what it means.

A degree-$d$ polynomial is fixed by its values at any $d+1$ distinct points, and the sum-of-squares condition on it can be
imposed at any $d+1$ distinct points with any positive weights. So the second job can be done at other points without
changing the polynomial program. With `--encodingPoints`, `pmp2sdp` keeps the approximation points for job 1 and writes
the SDP at a second set of $d+1$ **encoding points** $y_j$ with scalings $w_j$:

- rows of $B$ and $c$: $w_j \tilde R(y_j)$, where $\tilde R$ is the interpolant of the approximation-point values,
  evaluated from those values by barycentric Lagrange interpolation (one $(d+1)\times(d+1)$ matrix per block, applied
  to every matrix entry), not from monomial coefficients;
- bilinear bases: the polynomials orthonormal for $\sum_j w_j \delta(x - y_j)$ and $\sum_j w_j y_j \delta(x - y_j)$,
  times $\sqrt{w_j}$ and $\sqrt{w_j y_j}$, made by the three-term (Stieltjes) recurrence on the sampled vectors, with
  full reorthogonalization, instead of a Cholesky factorization of monomial moments (which loses most of the
  precision when the points span many decades);
- the objective and normalization are those of the default encoding, to working precision: the whole conversion runs at
  the extra precision and is then truncated, so the objective (and blocks of degree 0) can
  differ from the default's output in the last bits; `normalization.json` is identical.

In exact arithmetic the feasible set and the optimum are the same as with the default encoding (with or without a
reduced prefactor: both encode the same interpolant $R$). What changes is the SDP's coordinates, and with them the
number of iterations SDPB needs and the precision it needs to reach a verdict.

`pmp_info.json` lists the encoding points and scalings as each block's sample points and scalings (they are the
points at which the SDP is written).

Blocks of degree 0 (one sample point, e.g. an isolated constant constraint) have nothing to re-encode and keep the
default encoding, whatever the rule (with `file`, their entries in the file are not used).

## The rules for the encoding points

For a block with $m$ kept (reduced-prefactor) poles $p_i \le 0$, $d+1$ points, and a scale $a > 0$:

- `--encodingPoints functionalMatched`: with $g(x) = \prod_i (x + |p_i|)\, (x + a)^{d+1-m}$, the points solve
  $N_g(y_k) = k - \tfrac12$ ($k = 1, \dots, d+1$), where
  $N_g(x) = \tfrac{2}{\pi}\big[\sum_i \arctan\sqrt{x/|p_i|} + (d+1-m)\arctan\sqrt{x/a}\big]$, and the scalings are
  $w_k = 1/P_g(y_k)$ with $P_g(x) = \tfrac{g(x)}{\pi}\big[\sum_i \tfrac{\sqrt{|p_i|}}{x + |p_i|} + \tfrac{(d+1-m)\sqrt a}{x+a}\big]$.
  This is a quadrature rule exact for all polynomials of degree $\le d$ against $dx/(\sqrt{x}\, g(x))$ on
  $[0, \infty)$. Poles at $x = 0$ are counted in the $(x+a)$ factor.
- `--encodingPoints chebyshev`: the Chebyshev points of $t = x/(x+a)$, $y_k = a \tan^2\frac{(2k-1)\pi}{4(d+1)}$, with
  $w_k = 1 / \big(\prod_i (y_k - p_i)\, (y_k + a)^{d-m}\big)$.
- `--encodingPoints file --encodingPointsFile <f>`: points and scalings for every block, from a JSON file: a PMP-like
  object `{"PositiveMatrixWithPrefactorArray": [{"samplePoints": [...], "sampleScalings": [...]}, ...]}` (other fields
  are ignored; e.g. another run's `pmp_info.json`-like output) or the array alone, by block index.
- `--encodingPoints prefactorMatched` (the default): the SDP is written at the approximation points, as without the
  option.

The scale: $a = A \cdot \operatorname{median}_i |p_i|$ with $A$ = `--encodingPointsScale` (default 1), or $10A$ for a block
without poles; or one $a$ per block from `--encodingScaleFile`, a JSON array indexed by block.

Both rules put most of the points far beyond the approximation points, on the half-line $x \ge 0$ where the positivity
condition holds.

## Precision

The step from the approximation points to the far encoding points cancels digits. `pmp2sdp` therefore runs the whole
conversion at `--precision` + `--encodingExtraPrecision` bits (default 384) and writes the SDP at `--precision`
(binary block files carry `--precision`, rounded up to whole GMP limbs as SDPB reports it, so that SDPB reads them at its
own precision). It prints an estimate of the largest loss,
$\max_j \log_{10}\big(\max_n \sum_k |l_k(y_j) V_{kn}| \,/\, \max_n |\sum_k l_k(y_j) V_{kn}|\big)$ over rows ($l_k$ the
Lagrange polynomials of the approximation points, $V$ the unscaled rows), and stops with an error if fewer than 32 of the
extra bits would be left. On 3d Ising programs (below) the printed loss was 45 to 71 decimal digits; 384 extra bits are
115 digits. Raise `--encodingExtraPrecision` if `pmp2sdp` stops.

The conversion costs more than the default: 1.6 to 2.7 times the default's `pmp2sdp` time in the runs below (20 to 200
seconds more on one 128-core node), small against the SDPB calls it shortens.

## When to use it

Try it for programs with many poles and a large degree per block, where SDPB needs a high precision or many iterations:
for example the 3d Ising bootstrap with the stress tensor at nmax = 14, where the default needs 832 bits. Compare with the
default encoding on a point whose verdict you know before relying on it, and choose the precision by such a comparison:
the encoding changes the number of digits SDPB loses per iteration, so the precision a problem needs must be measured,
not assumed.

Limitations:
- `approx_objective` has not been tried on a two-set SDP. `spectrum` has been checked on one small problem, SDPB's
  nmax = 6 end-to-end test (a primal-dual optimal problem at 768 bits, duality gap 1e-30). With the two-set encoding,
  the objectives agree with the default's to 30 digits, y to 6e-23, the spectrum's zeros to 2e-22 and its OPE
  coefficients to 2e-14 (relative). The default encoding itself, started on another central path, differs from its own
  first run by 2e-23, 5e-23 and 5e-15. The blocks' `error` estimates in `spectrum.json` are residuals of the
  encoded SDP and differ between encodings (1e-28 to 1e-33 in the default encoding, down to 1e-51 in the two-set one).
- The two-set encoding has not been combined with nodes at infinity (`--reducedPrefactorNodesAtInfinity`, a separate
  pull request that moves some interpolation conditions of reduced blocks to x = ∞). If both are merged, `pmp2sdp`
  must refuse the combination with a clear error until the second step handles the moved rows; the lab's development
  branch does so.
- `functionalMatched` needs at least as many points as a block has nonzero kept poles. A block with more (for example
  11 poles on 5 points, as in SDPB's `1d-duplicate-poles` test, or a 12-pole block of degree 10, also when
  `--maxNumPoles 4` leaves 3 points) stops `pmp2sdp` with the assertion "more nonzero kept poles than points". Use
  `chebyshev` (or reduce `--maxNumPoles` below the number of points) for such blocks.

## Measurements

Ising with the stress tensor, nmax = 14, one excluded point (a feasibility problem, 300 blocks, a 35-pole reduced
prefactor, 3,407 free variables), `--encodingPoints functionalMatched --encodingPointsScale 1`, SDPB on two 128-core
Expanse nodes, the same SDPB build for every arm. "Checked" below means that an independent re-analysis of the
recorded runs reproduced the number; "measured once" means one run per arm, not repeated.

| encoding, precision | iterations to the exclusion | wall time of the SDPB call | |
|---|---|---|---|
| default, 832 bits | 106 | 5,487 s | checked |
| two-set, 832 bits | 60 | 3,114 s | iterations checked, time measured once |
| two-set, 320 bits | 60 | 1,716 s (two runs: 1,716.1 and 1,716.9 s) | checked |
| two-set, 256 bits | fails: the Cholesky factorization of SDPB's matrix Q fails at the first iteration | | measured once |
| default, 320 bits | 106 | 2,983 s | checked; marginal: its functional agrees with the 832-bit default's to 5.8 digits only |

- A whole call is **3.15 to 3.20 times faster** with the two-set encoding at 320 bits than with the default at 832
  (checked by an independent re-analysis; the matched factors use whole default times inferred from per-iteration
  timings, and the 3.197 end is a cross-node comparison that includes a restart). About 1.77× comes from fewer iterations (106 against 60) and 1.83× from the lower precision.
- The default also reaches the exclusion at 320 bits, with less margin (5.8 digits of agreement against 9.7 for the
  two-set encoding at 320: 9.66 in the objective, a median of 9.80 in μ). Against the default at 320, the two-set encoding is **1.74×** faster a call (checked; one
  default-320 run).
- At a fixed precision the two encodings cost the same per iteration.
- With a scale per block (`--encodingScaleFile`) derived from the functional of a neighbouring point, the same point
  was excluded at 256 bits in 55 iterations (checked: an independent repeat reproduced the 55 iterations and the 7.58 digits of agreement).
- The number of iterations depends on the scale $a$: with $A = 0.5$ instead of 1, both rules needed 102 to 104
  iterations (measured once). The rule (`functionalMatched` or `chebyshev`) mattered less than the scale in these
  tests. Smaller problems gained less: at nmax = 10 (152 blocks, `--maxNumPoles 20`), 62 against 70 iterations at every
  precision from 256 to 832 bits (measured once).

These numbers come from a single excluded point at one size. They are evidence for this kind of problem, not a
guarantee for others.
