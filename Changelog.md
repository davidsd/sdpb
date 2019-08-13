# Version 2.1.3

## sdpb

- Fixed a bug in load balancing that caused jobs to crash.

- Made some very, very small problems run faster.

## pvm2sdp

- Added the ability to read lists of input files from a file.

## sdp2input

- Added the ability to read lists of input files from a file.

- Significantly reduced memory use for runs with many cores.

# Version 2.1.2

## sdpb

- Significantly reduced memory use for large runs with many cores.

# Version 2.1.1

## sdpb

- Fixed a bug that was erroneously rejecting options in parameter
  files.

# Version 2.1.0

## sdpb

- Fixed the final solution output of the `x` vector, and enabled plain
  text output of the `X` and `Y` diagonal block matrices with the new
  `--writeMatrices` option.  This required splitting up these vectors
  and matrices into separate files.  So the option `outFile` has been
  changed to `outDir`, and all of the plain text solution files go
  into `outDir`.
  
- Added the ability to restart from plain text output generated with
  `--writeMatrices`.  This is useful if restarting using different
  precision or a different number of cores.  Due to limitations in
  GMP, a few bits of precision are lost during the translation to and
  from plain text.  So a run restarted from plain text will not be
  bitwise identical.

- Reorganized the options.

# Version 2.0.4

## sdpb

- Significantly reduced memory usage for large parallel jobs.  In
  addition to ordinary improvements in memory use, there is a new
  option: procGranularity.  This option can further reduce memory use at
  the cost of slower performance.

- Reordered how stopping conditions are checked.  The most visible
  change is that SDPB now checks for primal solutions before primal
  jumps.  This also fixes a bug where SDPB would hang with the option
  --findPrimalFeasible.

- Changed the logic for how output and checkpoint directory names are
  generated.  Instead of replacing any extension on sdpDir, sdpb now
  unconditionally adds ".out" or ".ck".  This prevents problems when,
  for example, there is a decimal number in sdpDir.

- Added a --version option.

## pvm2sdp

- Significantly reduced memory usage.

## sdp2input

- Fixed a corner case that was rejecting valid polynomials.

# Version 2.0.3

## sdpb

- Fix computation of the primal error

# Version 2.0.2

## sdpb

- Fix a bug that was erroneously writing checkpoints during timing
  runs

## sdp2input

- Fix a bug when computing normalization
- Handle parsing of constant polynomials in Mathematica input
- Handle DampedRational's that are just a constant

# Version 2.0.1

## sdpb

- Include training time in maxRunttime
- Fix a corner case in boost autodetection

# Version 2.0.0

- Initial release.
