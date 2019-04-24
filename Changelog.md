# Version 2.0.4

## sdpb

- Significantly reduced memory usage for large parallel jobs.

- Reordered how stopping conditions are checked.  The most visible
  change is that SDPB now checks for primal solutions before primal
  jumps.  This also fixes a bug where SDPB would hang with the option
  --findPrimalFeasible.

- Changed the logic for how output and checkpoint directory names are
  generated.  Instead of replacing any extension on sdpDir, sdpb now
  unconditionally adds ".out" or ".ck".  This prevents problems when,
  for example, there is a decimal number in sdpDir.

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
