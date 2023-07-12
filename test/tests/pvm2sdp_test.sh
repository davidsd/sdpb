#!/bin/bash

# setup
source test/common_test_setup.sh || exit 1

# Usage:
# run_pvm2sdp input output
function run_pvm2sdp() {
  mpirun -n 2 ./build/pvm2sdp 1024 $@
}

pvm2sdp_out=$TEST_OUT_DIR/pvm2sdp/sdp.zip
TEST_RUN_SUCCESS "run pvm2sdp" run_pvm2sdp $TEST_DATA_DIR/pvm2sdp/file_list.nsv $pvm2sdp_out
TEST_RUN_SUCCESS "pvm2sdp check output" diff $pvm2sdp_out $TEST_DATA_DIR/sdp.zip

# create file and prohibit writing
function touch_no_write() {
  local filename=$1
  local dir_name=$(dirname "$filename")
  rm -rf "$dir_name"
  mkdir -p "$dir_name"
  touch "$filename"
  chmod a-w "$filename"
}

io_tests="$TEST_OUT_DIR/io_tests"

touch_no_write $io_tests/pvm2sdp/sdp.zip
TEST_RUN_FAILS "pvm2sdp cannot write zip" run_pvm2sdp $TEST_DATA_DIR/pvm2sdp/pvm.xml $io_tests/pvm2sdp/sdp.zip

file_does_not_exist="$io_tests/file_does_not_exist"
mkdir -p $file_does_not_exist
echo "file_does_not_exist" >$file_does_not_exist/file_list.nsv
TEST_RUN_FAILS "corrupt file_list.nsv" run_pvm2sdp $file_does_not_exist/file_list.nsv $file_does_not_exist/sdp.zip
