#!/bin/bash

# setup
source test/common_test_setup.sh || exit 1

input_dir=$TEST_DATA_DIR/spectrum
output_dir=$TEST_OUT_DIR/spectrum
rm -rf $output_dir

echo "================"
echo "Running spectrum tests..."
TEST_RUN_SUCCESS "run spectrum" mpirun -n 2 build/spectrum --input=$input_dir/pvm.xml --solution=$input_dir/solution --threshold=1e-10 --format=PVM --output=$output_dir/spectrum.json --precision=1024
TEST_RUN_SUCCESS "check spectrum output" diff $output_dir/spectrum.json $input_dir/spectrum_orig.json
