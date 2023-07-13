#!/bin/bash

# setup
source test/common_test_setup.sh || exit 1

input_dir=$TEST_DATA_DIR/sdp2input
output_dir=$TEST_OUT_DIR/sdp2input
rm -rf $output_dir

function run_sdp2input() {
  local filename=$1
  local args=${*:2}
  # note that input files contain 214 decimal digits (~710 binary)
  # to have same results for .m and json, we have to set precision <=710
  mpirun -n 2 build/sdp2input --precision=512 --input="$input_dir/$filename" --output="$output_dir/$filename/sdp.zip" $args
}

function sdp2input_run_test() {
  local filename=$1
  local args=${*:2}
  local result=$output_dir/$filename/sdp.zip
  local orig=$input_dir/sdp_orig.zip
  TEST_RUN_SUCCESS "run sdp2input $filename" run_sdp2input "$filename" $args
  TEST_RUN_SUCCESS "check sdp2input result for $filename" DIFF_ZIP_IGNORE_CONTROL "$result" "$orig"
}

sdp2input_run_test "sdp2input_test.json"
sdp2input_run_test "sdp2input_split.nsv"
sdp2input_run_test "sdp2input_test.m" --debug=true

TEST_RUN_SUCCESS "sdp2input profiling.0" CHECK_FILE_NOT_EMPTY "$output_dir/sdp2input_test.m/sdp.zip.profiling.0"
TEST_RUN_SUCCESS "sdp2input profiling.1" CHECK_FILE_NOT_EMPTY "$output_dir/sdp2input_test.m/sdp.zip.profiling.1"
