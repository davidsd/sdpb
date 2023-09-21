#!/bin/sh

if [[ $1 == "--help" ]]; then
  echo "Pass custom mpirun command line as arguments to the testing script, e.g.:"
  echo "  ./test/run_all_tests.sh mpirun"
  echo "  ./test/run_all_tests.sh srun --mpi=pmi2"
  echo "This command will be stored in MPI_RUN_COMMAND variable."
  echo "By default, MPI_RUN_COMMAND=mpirun"
  exit
fi

# Set working directory to sdpb root
script_path="$0" # sdpb/test/run_all_test.sh
cd "$(dirname "$script_path")" # sdpb/test/
cd .. # sdpb/
echo "root directory: $PWD"

# setup custom mpirun command, e.g. "srun --mpi=pmi2"
if [ $# -ge 1 ]; then
  MPI_RUN_COMMAND="$@"
  echo "MPI_RUN_COMMAND=$MPI_RUN_COMMAND"
else
  MPI_RUN_COMMAND=mpirun
fi

# Run unit_tests and integration_tests with timing info and custom mpirun command
# For more command-line options, see
# https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md

# unit_tests: run 2 processes
echo time $MPI_RUN_COMMAND -n 2 ./build/unit_tests --durations yes
time $MPI_RUN_COMMAND -n 2 ./build/unit_tests --durations yes || { exit $?; }

# integration_tests
echo time ./build/integration_tests --durations yes --mpirun="$MPI_RUN_COMMAND"
time ./build/integration_tests --durations yes --mpirun="$MPI_RUN_COMMAND" || { exit $?; }

echo "$0: ALL TESTS PASSED"