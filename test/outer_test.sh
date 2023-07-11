#!/bin/bash

# TODO merge with run_test.sh

# setup
common_test_setup=test/common_test_setup.sh
if [ ! -f $common_test_setup ];then
  echo "$common_test_setup not found. Run this script from sdpb root directory"
  exit 1
fi
source $common_test_setup

data_dir=$TEST_DATA_DIR/outer_limits
output_dir=$TEST_OUT_DIR/outer_limits
rm -rf $output_dir
mkdir -p $output_dir # TODO create dir in Outer_Parameters.cxx instead

TEST_RUN_SUCCESS "run outer_limits" ./build/outer_limits --functions $data_dir/toy_functions.json --out $output_dir/toy_functions_out.json --checkpointDir $output_dir/ck --points $data_dir/toy_functions_points.json --precision=128 --dualityGapThreshold=1e-10 --primalErrorThreshold=1e-10 --dualErrorThreshold=1e-10 --initialMatrixScalePrimal=1e1 --initialMatrixScaleDual=1e1 --maxIterations=1000 --verbosity=0
TEST_RUN_SUCCESS "check outer_limits output" diff $output_dir/toy_functions_out.json $data_dir/toy_functions_out_orig.json
