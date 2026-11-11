#!/bin/sh

usage() {
  echo "Usage: ./test/run_SDPLIB.sh [--build-dir=DIR] [MPI_RUN_COMMAND...]"
  echo "DIR contains integration_tests. By default, DIR=build in the sdpb root directory."
  echo "Pass custom mpirun command line as arguments to the testing script, e.g.:"
  echo "  ./test/run_SDPLIB.sh mpirun"
  echo "  ./test/run_SDPLIB.sh srun --mpi=pmi2"
  echo "This command will be stored in MPI_RUN_COMMAND variable."
  echo "By default, MPI_RUN_COMMAND=mpirun --oversubscribe"
}

# Options come before MPI_RUN_COMMAND.
# Resolve DIR before changing the working directory.
BUILD_DIR=""
while [ $# -ge 1 ]; do
  case "$1" in
    --help)
      usage
      exit
      ;;
    --build-dir=*)
      BUILD_DIR="$(cd "${1#--build-dir=}" && pwd)" || exit 1
      shift
      ;;
    --build-dir)
      [ $# -ge 2 ] || { usage; exit 1; }
      BUILD_DIR="$(cd "$2" && pwd)" || exit 1
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

# Set working directory to sdpb root
script_path="$0" # sdpb/test/run_SDPLIB.sh
cd "$(dirname "$script_path")" || exit 1 # sdpb/test/
cd .. # sdpb/
echo "root directory: $PWD"
[ -n "$BUILD_DIR" ] || BUILD_DIR="$PWD/build"
echo "build directory: $BUILD_DIR"

# setup custom mpirun command, e.g. "srun --mpi=pmi2"
if [ $# -ge 1 ]; then
  MPI_RUN_COMMAND="$*"
  echo "MPI_RUN_COMMAND=$MPI_RUN_COMMAND"
else
  MPI_RUN_COMMAND="mpirun --oversubscribe"
fi

# Run unit_tests and integration_tests with timing info and custom mpirun command
# For more command-line options, see
# https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md

# integration_tests
echo time "$BUILD_DIR/integration_tests" SDPLIB --durations yes --mpirun="$MPI_RUN_COMMAND"
time "$BUILD_DIR/integration_tests" SDPLIB --durations yes --mpirun="$MPI_RUN_COMMAND" || { exit $?; }

echo "$0: ALL TESTS PASSED"