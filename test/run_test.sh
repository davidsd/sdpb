#!/bin/bash

# Run this from the top level directory
result=0

function test_run_with_exit_code() {
  local expected_exit_code=${1}
  local name=${2}
  local cmd=${*:3}
  if [ $expected_exit_code == 0 ]; then
    $cmd >/dev/null
  else
    $cmd >/dev/null 2>/dev/null
  fi
  local cmd_exit_code=$?
  if [ $cmd_exit_code == $expected_exit_code ]; then
    echo "PASS $name"
  else
    echo "FAIL $name: exit code $cmd_exit_code, expected $expected_exit_code"
    result=1
  fi
}

function test_run_success() {
  test_run_with_exit_code 0 "$@"
}

function test_run_fails() {
  test_run_with_exit_code 1 "$@"
}

rm -rf test/test/
test_run_success "pvm2sdp" ./build/pvm2sdp 1024 test/file_list.nsv test/test

rm -f test/test.out
test_run_success "SDPB" ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test --verbosity=0
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/profile_error.profiling.0
chmod a-w test/io_tests/profile_error.profiling.0
test_run_fails "write profile" mpirun -n 2 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test -c test/io_tests/profile_error --verbosity=2 --maxIterations=1
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/control.json
chmod a-w test/io_tests/control.json
test_run_fails "control" mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/block_0.json
chmod a-w test/io_tests/block_0.json
test_run_fails "blocks write" mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests
rm -rf test/io_tests

mkdir -p test/io_tests
touch test/io_tests/objectives.json
chmod a-w test/io_tests/objectives.json
test_run_fails "objectives" mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/test.xml test/io_tests
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/out.txt
chmod a-w test/io_tests/out/out.txt
test_run_fails "out.txt" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/x_0.txt
chmod a-w test/io_tests/out/x_0.txt
test_run_fails "x_0.txt" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/y.txt
chmod a-w test/io_tests/out/y.txt
test_run_fails "y.txt" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/X_matrix_0.txt
chmod a-w test/io_tests/out/X_matrix_0.txt
test_run_fails "X_matrix_0.txt" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --writeSolution=x,y,X,Y --verbosity=0
rm -rf test/io_tests

mkdir -p test/io_tests
mkdir test/io_tests/out
touch test/io_tests/out/Y_matrix_0.txt
chmod a-w test/io_tests/out/Y_matrix_0.txt
test_run_fails "Y_matrix_0.txt" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --writeSolution=x,y,X,Y --verbosity=0
rm -rf test/io_tests

mkdir -p test/io_tests
echo "file_does_not_exist" >test/io_tests/file_list.nsv
test_run_fails "file_list.nsv" mpirun -n 1 --quiet ./build/pvm2sdp 1024 test/io_tests/file_list.nsv test/io_tests
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
perl -p -0777 -i -e 'substr($_,1138,1)^=chr(1<<5)' test/io_tests/test
test_run_fails "input corruption" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 #
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
mpirun -n 1 ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
chmod a-r test/io_tests/out/X_matrix_0.txt
test_run_fails "text checkpoint read" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/out -o test/io_tests/out_new --maxIterations=1 --verbosity=0 #
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
mpirun -n 1 ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
rm test/io_tests/out/X_matrix_0.txt
touch test/io_tests/out/X_matrix_0.txt
test_run_fails "text checkpoint header" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/out -o test/io_tests/out_new --maxIterations=1 --verbosity=0
rm -rf test/io_tests

mkdir -p test/io_tests
cp -r test/test test/io_tests
mpirun -n 1 ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/ck -o test/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
head -n 2 test/io_tests/out/X_matrix_0.txt >test/io_tests/out/X_matrix_0.txt
test_run_fails "text checkpoint data" mpirun -n 1 --quiet ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s test/io_tests/test -c test/io_tests/out -o test/io_tests/out_new --maxIterations=1 --verbosity=0
rm -rf test/io_tests

exit $result
