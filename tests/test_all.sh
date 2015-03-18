#!/bin/bash

# Tobias Wood 2015
# Driver for all tests in order

source ./test_common.sh

run_test ./test_nifti.sh
run_test ./test_despot.sh

