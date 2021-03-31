#!/bin/bash

# Run this from the top level directory
result=0

rm -rf test/test/
./build/pvm2sdp 1024 test/file_list.nsv test/test/
if [ $? == 0 ]
then
    echo "PASS pvm2sdp"
else
    echo "FAIL pvm2sdp"
    result=1
fi

rm -f test/test.out
./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ --verbosity=0
diff test/test_out test/test_out_orig
if [ $? == 0 ]
then
    echo "PASS SDPB"
else
    echo "FAIL SDPB"
    result=1
fi
rm -rf test/io_tests

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
touch test/io_tests/control.json
chmod a-w test/io_tests/control.json
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS control"
else
    echo "FAIL control"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/block_0.json
chmod a-w test/io_tests/block_0.json
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS blocks write"
else
    echo "FAIL blocks write"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/objectives.json
chmod a-w test/io_tests/objectives.json
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS objectives"
else
    echo "FAIL objectives"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/out.txt
chmod a-w test/io_tests/out/out.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
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
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
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
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
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
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --writeSolution=x,y,X,Y --verbosity=0 2>/dev/null
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
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test/ -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --writeSolution=x,y,X,Y --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS Y_matrix_0.txt"
else
    echo "FAIL Y_matrix_0.txt"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
echo "file_does_not_exist" > test/io_tests/file_list.nsv
mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/io_tests/file_list.nsv test/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS file_list.nsv"
else
    echo "FAIL file_list.nsv"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
chmod a-r test/io_tests/test/block_0.json
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS blocks read"
else
    echo "FAIL blocks read"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
rm test/io_tests/test/block_0.json
touch test/io_tests/test/block_0.json
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS blocks corruption"
else
    echo "FAIL blocks corruption"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
chmod a-r test/io_tests/test/control.json
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS control read"
else
    echo "FAIL control read"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
rm test/io_tests/test/control.json
touch test/io_tests/test/control.json
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS control empty"
else
    echo "FAIL control empty"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
chmod a-r test/io_tests/test/objectives.json
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS objectives read"
else
    echo "FAIL objectives read"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
rm test/io_tests/test/objectives.json
touch test/io_tests/test/objectives.json
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS objectives corrupted"
else
    echo "FAIL objectives corrupted"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
mpirun -n 1 ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
chmod a-r test/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/out -o test/io_tests/out_new --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS text checkpoint read"
else
    echo "FAIL text checkpoint read"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
mpirun -n 1 ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
rm test/io_tests/out/X_matrix_0.txt
touch test/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/out -o test/io_tests/out_new --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS text checkpoint header"
else
    echo "FAIL text checkpoint header"
    result=1
fi
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
mpirun -n 1 ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
head -n 2 test/io_tests/out/X_matrix_0.txt > test/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/out -o test/io_tests/out_new --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS text checkpoint data"
else
    echo "FAIL text checkpoint data"
    result=1
fi
rm -rf test/io_tests

exit $result
