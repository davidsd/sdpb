#!/bin/sh

if [ "$1" = "--help" ]; then
  echo "Pass custom mpirun command line as arguments to the testing script, e.g.:"
  echo "  ./test/run_SDPLIB.sh mpirun"
  echo "  ./test/run_SDPLIB.sh srun --mpi=pmi2"
  echo "This command will be stored in MPI_RUN_COMMAND variable."
  echo "By default, MPI_RUN_COMMAND=mpirun --oversubscribe"
  exit
fi

# Set working directory to sdpb root
script_path="$0" # sdpb/test/run_all_test.sh
cd "$(dirname "$script_path")" || exit 1 # sdpb/test/
cd .. # sdpb/
echo "root directory: $PWD"

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
echo time ./build/integration_tests SDPLIB --durations yes --mpirun="$MPI_RUN_COMMAND"
time ./build/integration_tests SDPLIB --durations yes --mpirun="$MPI_RUN_COMMAND" || { exit $?; }

echo "$0: ALL TESTS PASSED"