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
    echo "PASS Profile Error"
else
    echo "FAIL Profile Error"
    result=1
fi
rm -f test/io_tests/profile_error.profiling.0

exit $result
