#!/bin/bash

# Run this from the top level directory

rm -f test/test.out
./build/sdpb --choleskyStabilizeThreshold=1e-1000 --precision=1024 --maxThreads=4 --noFinalCheckpoint -s test/test.xml > /dev/null
diff test/test.out test/test.out.orig
if [ $? == 0 ]
then
    echo "PASS"
else
    echo "FAIL"
fi
    
