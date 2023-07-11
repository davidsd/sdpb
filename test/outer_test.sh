#!/bin/bash

# TODO merge with run_test.sh

data_dir=test/data/outer_limits

output_dir=test/out/outer_limits
rm -rf $output_dir
mkdir -p $output_dir # TODO create dir in Outer_Parameters.cxx instead

./build/outer_limits --functions $data_dir/toy_functions.json --out $output_dir/toy_functions_out.json --checkpointDir $output_dir/ck --points $data_dir/toy_functions_points.json --precision=128 --dualityGapThreshold=1e-10 --primalErrorThreshold=1e-10 --dualErrorThreshold=1e-10 --initialMatrixScalePrimal=1e1 --initialMatrixScaleDual=1e1 --maxIterations=1000 --verbosity=0
diff $output_dir/toy_functions_out.json $data_dir/toy_functions_out_orig.json
if [ $? == 0 ]; then
    echo "PASS outer_limits"
else
    echo "FAIL outer_limits"
    exit 1
fi
