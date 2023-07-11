#!/bin/bash

# TODO merge with run_test.sh

data_dir=test/data/spectrum

output_dir=test/out/spectrum
rm -rf $output_dir

# NB: solution is copied from test/out/sdpb/out
# Still doesn't work, see also https://github.com/davidsd/sdpb/issues/89
./build/spectrum --precision=128 -i $data_dir/spectrum_test.json --solution $data_dir/solution -o $output_dir/spectrum_test_spectrum.json --threshold=1e-20 --format=PMP
# TODO after fixing, check also diff
#diff $output_dir/spectrum_test_spectrum.json $data_dir/spectrum_test_spectrum_orig.json
if [ $? == 0 ]; then
    echo "PASS spectrum"
else
    echo "FAIL spectrum"
    exit 1
fi
