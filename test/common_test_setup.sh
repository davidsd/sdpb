#!/bin/bash

# Common functions and variables to be used for specific tests
# Add the following line at the beginning of your test script (e.g. run_all_tests.sh):
# source test/common_test_setup.sh

# only include once
[ -n "$TEST_DATA_DIR" ] && return

TEST_RESULT=0
TEST_FAILED_LIST=""

echo "================"
echo "Common test setup..."

# prepare main directories
TEST_DATA_DIR=test/data
if [ ! -d "$TEST_DATA_DIR" ]; then
  echo "$TEST_DATA_DIR directory does not exist. Run this script from sdpb root directory"
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

  local log_stdout="$TEST_LOG_DIR/$name.stdout.log"
  local log_stderr="$TEST_LOG_DIR/$name.stderr.log"

  echo "test_run_with_exit_code $expected_exit_code $name $cmd" >>"$log_stdout"

  if [ $expected_exit_code == 0 ]; then
    $cmd >>"$log_stdout"
  else
    $cmd >>"$log_stdout" 2>>"$log_stderr"
  fi
  local cmd_exit_code=$?
  if [ $cmd_exit_code == $expected_exit_code ]; then
    echo "PASS $name"
  else
    echo "FAIL $name: exit code $cmd_exit_code, expected $expected_exit_code"
    echo "stdout: $log_stdout"
    echo "stderr: $log_stderr"
    TEST_FAILED_LIST="${TEST_FAILED_LIST} '$name'"
    TEST_RESULT=1
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

echo "================"
