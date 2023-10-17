# Version 2.6.0

- New **INCOMPATIBLE** format for sdp.zip.
  Instead of a single `block_XXX.json` file now we use `block_info_XXX.json` and `block_data_XXX.json` (
  or `block_data_XXX.bin`), see [SDPB_input_format.md](docs/SDPB_input_format.md).
  See [#114](https://github.com/davidsd/sdpb/pull/114).

- Compact binary format instead of JSON for `block_data_XXX` in sdp.zip.
  When generating sdp.zip by `sdp2input` or `pvm2sdp`, you can use optional command-line argument to choose between
  binary and JSON formats:

```
sdp2input --outputFormat FORMAT
pvm2sdp FORMAT PRECISION INPUT... OUTPUT
```

where `FORMAT` is either `bin` (used by default) or `json`.
We recommend using the new binary format since it is more compact and efficient.
See [#114](https://github.com/davidsd/sdpb/pull/114), [#128](https://github.com/davidsd/sdpb/pull/128), [#119](https://github.com/davidsd/sdpb/pull/119), [#149](https://github.com/davidsd/sdpb/pull/149).

- New lightweight Docker image does not contain `scalar_blocks` and `blocks_3d` anymore.
  Separate images for
  [sdpb](https://hub.docker.com/r/bootstrapcollaboration/sdpb/tags),
  [scalar_blocks](https://hub.docker.com/r/bootstrapcollaboration/scalar_blocks/tags)
  and [blocks_3d](https://hub.docker.com/r/bootstrapcollaboration/blocks_3d/tags)
  are now uploaded to https://hub.docker.com/u/bootstrapcollaboration.
  See [#130](https://github.com/davidsd/sdpb/pull/130).

- CI/CD pipelines on [CircleCI](https://app.circleci.com/pipelines/github/davidsd/sdpb).
  Tests are run for each `git push` to the repo and for each pull request.
  Docker images for the `master` branch and for each release starting from 2.6.0
  are uploaded to [Docker Hub](https://hub.docker.com/r/bootstrapcollaboration/sdpb/tags) automatically.
  See [#133](https://github.com/davidsd/sdpb/pull/133),
  [#136](https://github.com/davidsd/sdpb/pull/136).

- Tests are reorganized and rewritten from shell scripts to C++ [Catch2](https://github.com/catchorg/Catch2) framework.
  Added unit tests and realistic end-to-end tests.
  See [#91](https://github.com/davidsd/sdpb/pull/91),
  [#102](https://github.com/davidsd/sdpb/pull/102),
  [#109](https://github.com/davidsd/sdpb/pull/109),
  [#119](https://github.com/davidsd/sdpb/pull/119).

- Switched from C++14 to C++17.
  See [#118](https://github.com/davidsd/sdpb/pull/118),
  [#121](https://github.com/davidsd/sdpb/pull/121).

- Updated installations
  for [BU](docs/site_installs/Boston.md), [Caltech](docs/site_installs/Caltech.md), [Expanse](docs/site_installs/Expanse.md),
  and [Harvard](docs/site_installs/Harvard.md) clusters.

See https://github.com/davidsd/sdpb/releases/tag/2.6.0 for the full changelog.

# Version 2.5.1

## outer_limits

- Treat the behavior near zero with special care.  This means that
  input files now need to specify `epsilon_value` in addition to
  `infinity_value`.

- Added the `--meshThreshold` option to customize how finely the
  mesh is divided when searching for negative regions.

- Added the `--useSVD` option to control whether to regularize the
  problem with an SVD.

- Added the `sdp2functions` and `pvm2functions` programs to convert
  JSON and XML input files to the format that `outer_limits` expects.

- Fix an uninitialized matrix bug.

- Increase the maximum number of allowed SVD iterations because we are
  working at high precision.

## spectrum

- Initial release.

# Version 2.5.0

## sdp2input & pvm2sdp

- Uses a new, **INCOMPATIBLE** JSON-based format for writing output.
  The output is written as a single zip file, with entries that are
  each zip files.  This reduces problems with running out of inodes on
  some filesystems.  The zip format also transparently catches file
  corruption.  The format is documented in
  [SDPB_input_format.md](docs/SDPB_input_format.md).

- Mathematica polynomial input can now handle a plain 'x' as a
  polynomial term, rather than requiring '1.0*x'.

## sdp2input

- Fixed a bug causing incorrect results when poles are closer together
  than 1e-2.

- Significantly reduced memory use when using multiple cores

## sdpb

- Uses the new incompatible json-based format from pvm2sdp and sdp2input.

- Can read input as a directory or a wide variety of archive formats
  (e.g. zip, tar, 7z).  This requires a new dependency on the
  libarchive library.

- Added the options --minPrimalStep and --minDualStep as stopping conditions.

- Clean up output by only printing the Block Grid Mapping when verbosity=2.

## outer_limits

- Initial release.

## approx_objective

- Initial release.

# Version 2.4.0

## sdpb

- Made checkpointing more robust.  This required a modification to the
  checkpoint format.  Old checkpoints can still be read in, but new
  checkpoints will be written out in the new format.  Old checkpoints
  will not be automatically deleted.

- Better handle errors when reading blocks.  SDPB now gives
  informative error messages instead of inscrutable MPI crash reports.

## sdp2input

- Added support for JSON input and removed support for XML input.
  Mathematica input is still supported.  As a result,
  [RapidJSON](http://rapidjson.org/) is now a build requirement.

# Version 2.3.1

## sdpb

- Fixed a race condition when performing a timing run that caused
  incorrect answers.

- Fixed leaks of MPI communicators SDPB and Elemental that caused
  crashes when running long jobs with Intel's MPI.

## sdp2input

- Fixed a corner case when reading Mathematica files that do not have
  precision markers.

# Version 2.3.0

## sdpb

- Replaced the option --writeMatrices with a more general option
  --writeSolution.  This enables fine grained control over which output
  files to write.

- Changed the default for --maxRuntime to effectively infinite (2^63
  seconds)

- Use matrix size as a heuristic when allocating cores for a timing
  run.  This significantly decreases memory used during the timing
  run, making it more comparable to the regular solver run.

- Write fewer files when there are constant constraints.

- Fix a corner case when running a small problem on many cores

## pvm2sdp

- Implemented a better check for bad elements when reading in data.
  This also happened to speed it up by about 15%.

- No longer needlessly promote constant constraints to linear
  constraints.

# Version 2.2.0

## sdpb

- Fixed some race conditions when writing files and running on large
  numbers of cores.

- Added a number of checks for errors when reading and writing files.
  This required modifying the text checkpoints format to include a
  newline at the end, so text checkpoint written with older versions
  of SDPB will not load with this new version.  Binary checkpoints are
  unaffected and will continue to work.

- Require a valid checkpoint if the initialCheckpointDir option is set.

# Version 2.1.3

## sdpb

- Fixed a bug in load balancing that caused jobs to crash on some
  machines.

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
