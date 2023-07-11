#!/bin/bash

# TODO merge with run_test.sh

# setup
common_test_setup=test/common_test_setup.sh
if [ ! -f $common_test_setup ];then
  echo "$common_test_setup not found. Run this script from sdpb root directory"
  exit 1
fi
source $common_test_setup

data_dir=$TEST_DATA_DIR/spectrum
output_dir=$TEST_OUT_DIR/spectrum
rm -rf $output_dir

# NB: solution is copied from test/out/sdpb/out
# Still doesn't work, see also https://github.com/davidsd/sdpb/issues/89
TEST_RUN_SUCCESS "run spectrum" ./build/spectrum --precision=128 -i $data_dir/spectrum_test.json --solution $data_dir/solution -o $output_dir/spectrum_test_spectrum.json --threshold=1e-20 --format=PMP
# TODO after fixing, check also diff
 TEST_RUN_SUCCESS "check spectrum result" diff $output_dir/spectrum_test_spectrum.json $data_dir/spectrum_test_spectrum_orig.json
