#!/bin/bash

# Common functions and variables to be used for specific tests
# Add the following line at the beginning of your test script (e.g. run_all_tests.sh):
# source test/common_test_setup.sh

# only include once
[ -n "$TEST_DATA_DIR" ] && return

TEST_SCRIPT_PATH="$0"

TEST_FAILED_LIST=""
TEST_PASSED_COUNT=0
TEST_FAILED_COUNT=0

# usage:
# echo -e "${SET_COLOR_RED}Red text${UNSET_COLOR}"
SET_COLOR_RED='\033[0;31m'
SET_COLOR_GREEN='\033[0;32m'
UNSET_COLOR='\e[0m'

function ECHO_RED() {
  local text="$*"
  echo -e "${SET_COLOR_RED}${text}${UNSET_COLOR}"
}
function ECHO_GREEN() {
  local text="$*"
  echo -e "${SET_COLOR_GREEN}${text}${UNSET_COLOR}"
}

echo "================"
echo "Common test setup..."

if [[ $1 == "--help" ]]; then
  echo "common_test_setup.sh help:"
  echo "Pass custom mpirun command line as arguments to the testing script, e.g.:"
  echo "  ./run_all_tests.sh mpirun -n 2"
  echo "  ./run_all_tests.sh srun -n \$SLURM_NTASKS --mpi=pmi2"
  echo "This command will be stored in MPI_RUN_COMMAND variable."
  echo "By default, MPI_RUN_COMMAND=mpirun -n 2"
  echo "Note that some tests require at least two MPI processes."
  exit
fi

# setup custom mpirun command, e.g. "srun -n $SLURM_NTASKS --mpi=pmi2"
if [ $# -ge 1 ]; then
  MPI_RUN_COMMAND="$@"
else
  MPI_RUN_COMMAND="mpirun -n 2"
fi

$MPI_RUN_COMMAND echo "test MPI_RUN_COMMAND" >/dev/null || {
  ECHO_RED "Invalid MPI_RUN_COMMAND=$MPI_RUN_COMMAND"
  exit 1
}
echo "MPI_RUN_COMMAND=$MPI_RUN_COMMAND"

# prepare main directories
TEST_DATA_DIR=test/data
if [ ! -d "$TEST_DATA_DIR" ]; then
  ECHO_RED "$PWD/$TEST_DATA_DIR: directory does not exist."
  exit 1
fi
echo "TEST_DATA_DIR=$TEST_DATA_DIR"

TEST_OUT_DIR=test/out
echo "TEST_OUT_DIR=$TEST_OUT_DIR"
rm -rf $TEST_OUT_DIR
mkdir -p $TEST_OUT_DIR

TEST_LOG_DIR=$TEST_OUT_DIR/log
echo "TEST_LOG_DIR=$TEST_LOG_DIR"
rm -rf $TEST_LOG_DIR
mkdir -p $TEST_LOG_DIR

# test functions

function TEST_RUN_WITH_EXIT_CODE() {
  local expected_exit_code=${1}
  local name=${2}
  local cmd=${*:3}

  # get current test script name, e.g. sdpb_test.sh
  local curr_test_filename=$CURR_TEST_PATH
  [[ -z $curr_test_filename ]] && curr_test_filename=$TEST_SCRIPT_PATH
  curr_test_filename=$(basename $curr_test_filename)

  # log paths
  # To ensure valid log filename, keep only alphanumerics and '_-.':
  local name_sanitized=${name//[^A-Za-z0-9_-.]/_}
  local curr_log_dir="$TEST_LOG_DIR/$curr_test_filename"
  mkdir -p $curr_log_dir
  local path_sanitized="$curr_log_dir/$name_sanitized"
  local log_stdout="$path_sanitized.stdout.log"
  local log_stderr="$path_sanitized.stderr.log"

  # name for printing
  name="$curr_test_filename/$name"

  # run and print to logs
  echo "$cmd" >>"$log_stdout"
  $cmd >>"$log_stdout" 2>>"$log_stderr"
  local cmd_exit_code=$?
  cat "$log_stderr" >>"$log_stdout"

  if [ $cmd_exit_code == $expected_exit_code ]; then
    ECHO_GREEN "PASS $name"
    ((TEST_PASSED_COUNT++))
  else
    echo "$name"
    echo "$cmd"
    cat "$log_stderr"
    ECHO_RED "FAIL $name: exit code $cmd_exit_code, expected $expected_exit_code"
    echo "stdout: $log_stdout"
    echo "stderr: $log_stderr"
    echo "----------------"
    TEST_FAILED_LIST="${TEST_FAILED_LIST} '$name'"
    ((TEST_FAILED_COUNT++))
  fi
}

function TEST_RUN_SUCCESS() {
  TEST_RUN_WITH_EXIT_CODE 0 "$@"
}

# TODO actually we should check not exit code, but error message from stderr!
# Program may fail for some other reason, and we'll never notice that.
function TEST_RUN_FAILS() {
  TEST_RUN_WITH_EXIT_CODE 1 "$@"
}

function CHECK_FILE_NOT_EMPTY() {
  filename="$1"
  if [[ -s "$filename" ]]; then
    return 0
  else
    return 1
  fi
}

function ZIP_SUMMARY() {
  local filename="$1"
  (unzip -vqq "$filename" | awk '{$2=""; $3=""; $4=""; $5=""; $6=""; print}' | sort -k3 -f)
}

# compare two sdp.zip archives by content ignoring control.json (it contains command which can be different)
# explanation: https://stackoverflow.com/a/61113635/3270684
function DIFF_ZIP_IGNORE_CONTROL() {
  local x=$1
  local y=$2
  diff --exclude=control.json \
    <(ZIP_SUMMARY "$x" | grep -v "control.json") \
    <(ZIP_SUMMARY "$y" | grep -v "control.json")
}

echo "================"
