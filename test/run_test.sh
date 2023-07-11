#!/bin/bash

# Run this from the top level directory
result=0

# prepare main directories
data_dir=test/data

out_dir=test/out
rm -rf $out_dir
mkdir -p $out_dir

log_dir=$out_dir/log
rm -rf $log_dir
mkdir -p $log_dir

pvm2sdp_out=$out_dir/pvm2sdp/sdp.zip
sdp_path=$data_dir/sdp.zip

failed_tests=""

# test functions

function test_run_with_exit_code() {
  local expected_exit_code=${1}
  local name=${2}
  local cmd=${*:3}

  local log_stdout="$log_dir/$name.stdout.log"
  local log_stderr="$log_dir/$name.stderr.log"

  echo "test_run_with_exit_code $expected_exit_code $name $cmd" >>"$log_stdout"

  if [ $expected_exit_code == 0 ]; then
    $cmd >>"$log_stdout"
  else
    $cmd >>"$log_stdout" 2>>"$log_stderr"
  fi
  local cmd_exit_code=$?
  if [ $cmd_exit_code == $expected_exit_code ]; then
    echo "PASS $name"
  else
    echo "FAIL $name: exit code $cmd_exit_code, expected $expected_exit_code"
    echo "stdout: $log_stdout"
    echo "stderr: $log_stderr"
    failed_tests="${failed_tests} '$name'"
    result=1
  fi
}

function test_run_success() {
  test_run_with_exit_code 0 "$@"
}

# TODO actually we should check not exit code, but error message from stderr!
# Program may fail for some other reason, and we'll never notice that.
function test_run_fails() {
  test_run_with_exit_code 1 "$@"
}

# shortcuts for running pvm2sdp & sdpb with common options

# Usage:
# run_pvm2sdp input output
function run_pvm2sdp() {
  mpirun -n 2 ./build/pvm2sdp 1024 $@
}

# set default parameters except for paths
# Usage:
# run_sdpb -s "sdpPath" -c "checkpointDir" -o "outDir"
function run_sdpb() {
  mpirun -n 2 ./build/sdpb --precision=1024 --noFinalCheckpoint --procsPerNode=1 $@
}

# set default parameters except for paths
# Usage:
# run_sdpb_default_sdp_custom_output_prefix "output_prefix"
function run_sdpb_default_sdp_custom_output_prefix() {
  local output_prefix=$1
  run_sdpb -s $sdp_path -c "$output_prefix/ck" -o "$output_prefix/out" "${*:2}"
}

echo "================"
echo "Running tests..."

test_run_success "run pvm2sdp" run_pvm2sdp $data_dir/pvm2sdp/file_list.nsv $pvm2sdp_out
test_run_success "pvm2sdp check output" diff $pvm2sdp_out $sdp_path
test_run_success "run sdpb" run_sdpb_default_sdp_custom_output_prefix $out_dir/sdpb
test_run_success "sdpb check result" diff $out_dir/sdpb/out $data_dir/sdpb/test_out_orig

# create file and prohibit writing
function touch_no_write() {
  local filename=$1
  local dir_name=$(dirname "$filename")
  rm -rf "$dir_name"
  mkdir -p "$dir_name"
  touch "$filename"
  chmod a-w "$filename"
}

echo "================"
echo "Running IO failure tests..."

io_tests="$out_dir/io_tests"

touch_no_write $io_tests/write_profile/ck.profiling.0
test_run_fails "write profile" run_sdpb_default_sdp_custom_output_prefix "$io_tests/write_profile" --maxIterations=1 --verbosity=2

touch_no_write $io_tests/pvm2sdp/sdp.zip
test_run_fails "pvm2sdp cannot write zip" run_pvm2sdp $data_dir/pvm2sdp/test.xml $io_tests/pvm2sdp/sdp.zip

function test_sdpb_nowrite() {
  local name="$1"
  local output_prefix="$io_tests/$name"
  touch_no_write "$output_prefix/out/$name"
  test_run_fails "nowrite $name" run_sdpb_default_sdp_custom_output_prefix "$output_prefix" --maxIterations=1 --writeSolution=x,y,X,Y
}

test_sdpb_nowrite "out.txt"
test_sdpb_nowrite "x_0.txt"
test_sdpb_nowrite "y.txt"
test_sdpb_nowrite "X_matrix_0.txt"
test_sdpb_nowrite "Y_matrix_0.txt"

file_does_not_exist="$io_tests/file_does_not_exist"
mkdir -p $file_does_not_exist
echo "file_does_not_exist" >$file_does_not_exist/file_list.nsv
test_run_fails "corrupt file_list.nsv" run_pvm2sdp $file_does_not_exist/file_list.nsv $file_does_not_exist/sdp.zip

input_corruption="$io_tests/input_corruption"
mkdir -p "$input_corruption"
cp -r $sdp_path "$input_corruption"
perl -p -0777 -i -e 'substr($_,1138,1)^=chr(1<<5)' $input_corruption/sdp.zip
test_run_fails "corrupt sdp.zip" run_sdpb -s $input_corruption/sdp.zip -c $input_corruption/ck -o $input_corruption/out --maxIterations=1

checkpoint_read="$io_tests/checkpoint_read"
run_sdpb_default_sdp_custom_output_prefix "$checkpoint_read" --maxIterations=1 --writeSolution=x,y,X,Y >/dev/null
chmod a-r $checkpoint_read/out/X_matrix_0.txt
test_run_fails "noread checkpoint X_matrix_0.txt" run_sdpb -s $sdp_path -c $checkpoint_read/out -o $checkpoint_read/out_new --maxIterations=1

checkpoint_header=$io_tests/checkpoint_header
run_sdpb_default_sdp_custom_output_prefix "$checkpoint_header" --maxIterations=1 --writeSolution=x,y,X,Y >/dev/null
rm $checkpoint_header/out/X_matrix_0.txt
touch $checkpoint_header/out/X_matrix_0.txt
test_run_fails "empty text checkpoint X_matrix_0.txt" run_sdpb -s $sdp_path -c $checkpoint_header/out -o $checkpoint_header/out_new --maxIterations=1

checkpoint_data=$io_tests/checkpoint_data
run_sdpb_default_sdp_custom_output_prefix $checkpoint_data --maxIterations=1 --writeSolution=x,y,X,Y >/dev/null
head -n 2 $checkpoint_data/out/X_matrix_0.txt >$checkpoint_data/out/X_matrix_0.txt
test_run_fails "corrupt text checkpoint X_matrix_0.txt" run_sdpb -s $sdp_path -c $checkpoint_data/out -o $checkpoint_data/out_new --maxIterations=1

echo "================"
if [ $result != 0 ]; then
  echo "FAILED TESTS: $failed_tests"
else
  echo "PASSED ALL TESTS"
fi
exit $result
