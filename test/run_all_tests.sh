#!/bin/bash

# Set working directory to sdpb root
script_path="$0" # sdpb/test/run_all_test.sh
cd "$(dirname "$script_path")" # sdpb/test/
cd .. # sdpb/

# setup
source test/common_test_setup.sh || { echo "Run this script from sdpb root directory"; exit 1; }

source test/sdp2input_test.sh
source test/sdpb_test.sh
source test/outer_test.sh
source test/spectrum_test.sh

echo "================"
if [ $TEST_RESULT != 0 ]; then
  echo "FAILED TESTS: $TEST_FAILED_LIST"
else
  echo "PASSED ALL TESTS"
fi
exit $TEST_RESULT
