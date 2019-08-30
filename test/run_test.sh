#!/bin/bash

# Run this from the top level directory

./build/pvm2sdp 1024 test/file_list.nsv test/test/
rm -f test/test.out
./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ > /dev/null
diff test/test_out test/test_out_orig
result=0
if [ $? == 0 ]
then
    echo "PASS Smoke test"
else
    echo "FAIL Smoke test"
    result=1
fi
    
mkdir -p test/io_tests

touch test/io_tests/profile_error.profiling.0
chmod a-w test/io_tests/profile_error.profiling.0
mpirun -n 2 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/profile_error --verbosity=2 --maxIterations=1 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS write profile"
else
    echo "FAIL write profile"
    result=1
fi
rm -f test/io_tests/profile_error.profiling.0

touch test/io_tests/bilinear_bases.0
chmod a-w test/io_tests/bilinear_bases.0
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS bilinear_bases"
else
    echo "FAIL bilinear_bases"
    result=1
fi
rm -f test/io_tests/bilinear_bases.0

exit $result
