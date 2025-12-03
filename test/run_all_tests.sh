#!/bin/sh

if [ "$1" = "--help" ]; then
  echo "Pass custom mpirun command line as arguments to the testing script, e.g.:"
  echo "  ./test/run_all_tests.sh mpirun"
  echo "  ./test/run_all_tests.sh srun --mpi=pmi2"
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

# Find existing time command, see https://stackoverflow.com/a/677212/3270684
if command -v /usr/bin/time >/dev/null 2>&1; then
    TIME_CMD="/usr/bin/time"
elif command -v time >/dev/null 2>&1; then
    TIME_CMD="time"
else
    TIME_CMD=""
fi

# Run unit_tests and integration_tests with timing info and custom mpirun command
# For more command-line options, see
# https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md

# NB: for calculate_matrix_square_test to pass, number of processes must be a multiple of 6
echo $TIME_CMD $MPI_RUN_COMMAND -n 6 ./build/unit_tests --durations yes
$TIME_CMD $MPI_RUN_COMMAND -n 6 ./build/unit_tests --durations yes || { exit $?; }

# integration_tests
echo $TIME_CMD ./build/integration_tests --durations yes --mpirun="$MPI_RUN_COMMAND"
$TIME_CMD ./build/integration_tests --durations yes --mpirun="$MPI_RUN_COMMAND" || { exit $?; }

echo "$0: ALL TESTS PASSED"