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
rm -rf test/io_tests

mkdir -p test/io_tests
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
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/blocks.0
chmod a-w test/io_tests/blocks.0
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS blocks"
else
    echo "FAIL blocks"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/free_var_matrix.0
chmod a-w test/io_tests/free_var_matrix.0
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS free_var_matrix"
else
    echo "FAIL free_var_matrix"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/objectives
chmod a-w test/io_tests/objectives
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS objectives"
else
    echo "FAIL objectives"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/primal_objective_c.0
chmod a-w test/io_tests/primal_objective_c.0
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS primal_objective_c"
else
    echo "FAIL primal_objective_c"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/out.txt
chmod a-w test/io_tests/out/out.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS out.txt"
else
    echo "FAIL out.txt"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/x_0.txt
chmod a-w test/io_tests/out/x_0.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS x_0.txt"
else
    echo "FAIL x_0.txt"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/y.txt
chmod a-w test/io_tests/out/y.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS y.txt"
else
    echo "FAIL y.txt"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/X_matrix_0.txt
chmod a-w test/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --writeMatrices 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS X_matrix_0.txt"
else
    echo "FAIL X_matrix_0.txt"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/Y_matrix_0.txt
chmod a-w test/io_tests/out/Y_matrix_0.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --writeMatrices 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS Y_matrix_0.txt"
else
    echo "FAIL Y_matrix_0.txt"
    result=1
fi
rm -rf test/io_tests

exit $result
