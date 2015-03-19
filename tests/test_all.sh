#!/bin/bash

# Tobias Wood 2015
# Driver for all tests in order

source ./test_common.sh

run_test "NIFTI tests" ./test_nifti.sh
run_test "Single-Component DESPOT tests" ./test_scdespot.sh
run_test "Multi-Component DESPOT tests" ./test_mcdespot.sh

