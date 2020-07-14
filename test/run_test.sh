#!/bin/bash

# Run this from the top level directory
result=0

: "${SDPB_PVM2SDP:=./build/pvm2sdp}"
: "${SDPB_SDPB:=./build/sdpb}"
: "${SDPB_TESTPREFIX:=test}"

rm -rf ${SDPB_TESTPREFIX}/test/
${SDPB_PVM2SDP} 1024 ${SDPB_TESTPREFIX}/file_list.nsv ${SDPB_TESTPREFIX}/test/
if [ $? == 0 ]
then
    echo "PASS pvm2sdp"
else
    echo "FAIL pvm2sdp"
    result=1
fi

rm -f ${SDPB_TESTPREFIX}/test.out
${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/test/ --verbosity=0
diff ${SDPB_TESTPREFIX}/test_out ${SDPB_TESTPREFIX}/test_out_orig
if [ $? == 0 ]
then
    echo "PASS SDPB"
else
    echo "FAIL SDPB"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
touch ${SDPB_TESTPREFIX}/io_tests/profile_error.profiling.0
chmod a-w ${SDPB_TESTPREFIX}/io_tests/profile_error.profiling.0
mpirun -n 2 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/test/ -c ${SDPB_TESTPREFIX}/io_tests/profile_error --verbosity=2 --maxIterations=1 2>/dev/null > /dev/null
if [ $? != 0 ]
then
    echo "PASS write profile"
else
    echo "FAIL write profile"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
touch ${SDPB_TESTPREFIX}/io_tests/bilinear_bases.0
chmod a-w ${SDPB_TESTPREFIX}/io_tests/bilinear_bases.0
mpirun -n 1 --quiet ${SDPB_PVM2SDP} 1024 ${SDPB_TESTPREFIX}/test.xml ${SDPB_TESTPREFIX}/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS bilinear_bases"
else
    echo "FAIL bilinear_bases"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
touch ${SDPB_TESTPREFIX}/io_tests/blocks.0
chmod a-w ${SDPB_TESTPREFIX}/io_tests/blocks.0
mpirun -n 1 --quiet ${SDPB_PVM2SDP} 1024 ${SDPB_TESTPREFIX}/test.xml ${SDPB_TESTPREFIX}/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS blocks"
else
    echo "FAIL blocks"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
touch ${SDPB_TESTPREFIX}/io_tests/free_var_matrix.0
chmod a-w ${SDPB_TESTPREFIX}/io_tests/free_var_matrix.0
mpirun -n 1 --quiet ${SDPB_PVM2SDP} 1024 ${SDPB_TESTPREFIX}/test.xml ${SDPB_TESTPREFIX}/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS free_var_matrix"
else
    echo "FAIL free_var_matrix"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
touch ${SDPB_TESTPREFIX}/io_tests/objectives
chmod a-w ${SDPB_TESTPREFIX}/io_tests/objectives
mpirun -n 1 --quiet ${SDPB_PVM2SDP} 1024 ${SDPB_TESTPREFIX}/test.xml ${SDPB_TESTPREFIX}/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS objectives"
else
    echo "FAIL objectives"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
touch ${SDPB_TESTPREFIX}/io_tests/primal_objective_c.0
chmod a-w ${SDPB_TESTPREFIX}/io_tests/primal_objective_c.0
mpirun -n 1 --quiet ${SDPB_PVM2SDP} 1024 ${SDPB_TESTPREFIX}/test.xml ${SDPB_TESTPREFIX}/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS primal_objective_c"
else
    echo "FAIL primal_objective_c"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
mkdir ${SDPB_TESTPREFIX}/io_tests/out
touch ${SDPB_TESTPREFIX}/io_tests/out/out.txt
chmod a-w ${SDPB_TESTPREFIX}/io_tests/out/out.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/test/ -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS out.txt"
else
    echo "FAIL out.txt"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
mkdir ${SDPB_TESTPREFIX}/io_tests/out
touch ${SDPB_TESTPREFIX}/io_tests/out/x_0.txt
chmod a-w ${SDPB_TESTPREFIX}/io_tests/out/x_0.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/test/ -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS x_0.txt"
else
    echo "FAIL x_0.txt"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
mkdir ${SDPB_TESTPREFIX}/io_tests/out
touch ${SDPB_TESTPREFIX}/io_tests/out/y.txt
chmod a-w ${SDPB_TESTPREFIX}/io_tests/out/y.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/test/ -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS y.txt"
else
    echo "FAIL y.txt"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
mkdir ${SDPB_TESTPREFIX}/io_tests/out
touch ${SDPB_TESTPREFIX}/io_tests/out/X_matrix_0.txt
chmod a-w ${SDPB_TESTPREFIX}/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/test/ -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --writeSolution=x,y,X,Y --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS X_matrix_0.txt"
else
    echo "FAIL X_matrix_0.txt"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
mkdir ${SDPB_TESTPREFIX}/io_tests/out
touch ${SDPB_TESTPREFIX}/io_tests/out/Y_matrix_0.txt
chmod a-w ${SDPB_TESTPREFIX}/io_tests/out/Y_matrix_0.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/test/ -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --writeSolution=x,y,X,Y --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS Y_matrix_0.txt"
else
    echo "FAIL Y_matrix_0.txt"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
echo "file_does_not_exist" > ${SDPB_TESTPREFIX}/io_tests/file_list.nsv
mpirun -n 1 --quiet ${SDPB_PVM2SDP} 1024 ${SDPB_TESTPREFIX}/io_tests/file_list.nsv ${SDPB_TESTPREFIX}/io_tests 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS file_list.nsv"
else
    echo "FAIL file_list.nsv"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
chmod a-r ${SDPB_TESTPREFIX}/io_tests/test/blocks.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS blocks"
else
    echo "FAIL blocks"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
chmod a-r ${SDPB_TESTPREFIX}/io_tests/test/blocks.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS blocks"
else
    echo "FAIL blocks"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
rm ${SDPB_TESTPREFIX}/io_tests/test/blocks.0
touch ${SDPB_TESTPREFIX}/io_tests/test/blocks.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS blocks corruption"
else
    echo "FAIL blocks corruption"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
head -n 1 ${SDPB_TESTPREFIX}/io_tests/test/blocks.0 > ${SDPB_TESTPREFIX}/io_tests/test/blocks.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS vector size"
else
    echo "FAIL vector size"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
head -n 2 ${SDPB_TESTPREFIX}/io_tests/test/blocks.0 > ${SDPB_TESTPREFIX}/io_tests/test/blocks.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS vector element"
else
    echo "FAIL vector element"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
chmod a-r ${SDPB_TESTPREFIX}/io_tests/test/bilinear_bases.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS bilinear read"
else
    echo "FAIL bilinear read"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
rm ${SDPB_TESTPREFIX}/io_tests/test/bilinear_bases.0
touch ${SDPB_TESTPREFIX}/io_tests/test/bilinear_bases.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS bilinear empty"
else
    echo "FAIL bilinear empty"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
head -n 2 ${SDPB_TESTPREFIX}/io_tests/test/bilinear_bases.0 > ${SDPB_TESTPREFIX}/io_tests/test/bilinear_bases.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS bilinear header"
else
    echo "FAIL bilinear header"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
head -n 26 ${SDPB_TESTPREFIX}/io_tests/test/bilinear_bases.0 > ${SDPB_TESTPREFIX}/io_tests/test/bilinear_bases.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS bilinear data"
else
    echo "FAIL bilinear data"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
chmod a-r ${SDPB_TESTPREFIX}/io_tests/test/free_var_matrix.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS free_var_matrix read"
else
    echo "FAIL free_var_matrix read"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
rm ${SDPB_TESTPREFIX}/io_tests/test/free_var_matrix.0
touch ${SDPB_TESTPREFIX}/io_tests/test/free_var_matrix.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS free_var_matrix header"
else
    echo "FAIL free_var_matrix header"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
head -n 2 ${SDPB_TESTPREFIX}/io_tests/test/free_var_matrix.0 > ${SDPB_TESTPREFIX}/io_tests/test/free_var_matrix.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS free_var_matrix data"
else
    echo "FAIL free_var_matrix data"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
chmod a-r ${SDPB_TESTPREFIX}/io_tests/test/objectives
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS objectives read"
else
    echo "FAIL objectives read"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
rm ${SDPB_TESTPREFIX}/io_tests/test/objectives
touch ${SDPB_TESTPREFIX}/io_tests/test/objectives
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS objectives corrupted"
else
    echo "FAIL objectives corrupted"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
chmod a-r ${SDPB_TESTPREFIX}/io_tests/test/primal_objective_c.0
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS primal_objective_c"
else
    echo "FAIL primal_objective_c"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
mpirun -n 1 ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
chmod a-r ${SDPB_TESTPREFIX}/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/out -o ${SDPB_TESTPREFIX}/io_tests/out_new --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS text checkpoint read"
else
    echo "FAIL text checkpoint read"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
mpirun -n 1 ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
rm ${SDPB_TESTPREFIX}/io_tests/out/X_matrix_0.txt
touch ${SDPB_TESTPREFIX}/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/out -o ${SDPB_TESTPREFIX}/io_tests/out_new --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS text checkpoint header"
else
    echo "FAIL text checkpoint header"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

mkdir -p ${SDPB_TESTPREFIX}/io_tests
cp -r ${SDPB_TESTPREFIX}/test ${SDPB_TESTPREFIX}/io_tests
mpirun -n 1 ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/ck -o ${SDPB_TESTPREFIX}/io_tests/out --maxIterations=1 --verbosity=0 --writeSolution=x,y,X,Y
head -n 2 ${SDPB_TESTPREFIX}/io_tests/out/X_matrix_0.txt > ${SDPB_TESTPREFIX}/io_tests/out/X_matrix_0.txt
mpirun -n 1 --quiet ${SDPB_SDPB} --precision=1024 --noFinalCheckpoint --procsPerNode=1 -s ${SDPB_TESTPREFIX}/io_tests/test -c ${SDPB_TESTPREFIX}/io_tests/out -o ${SDPB_TESTPREFIX}/io_tests/out_new --maxIterations=1 --verbosity=0 2>/dev/null
if [ $? != 0 ]
then
    echo "PASS text checkpoint data"
else
    echo "FAIL text checkpoint data"
    result=1
fi
rm -rf ${SDPB_TESTPREFIX}/io_tests

exit $result
