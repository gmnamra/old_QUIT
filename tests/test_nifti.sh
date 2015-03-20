#!/bin/bash

# Tobias Wood 2015
# Simple test script for Nifti I/O tools

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="nifti_$QUITVER"
mkdir -p $DATADIR
cd $DATADIR

run_test "CREATE" $QUITDIR/niicreate blank.nii -d "16 16 16" -v "2 2 2"
run_test "READ" $QUITDIR/niihdr blank.nii
SILENCE_TEST="0"

cd ..
