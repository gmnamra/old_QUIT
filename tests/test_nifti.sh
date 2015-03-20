#!/bin/bash

# Tobias Wood 2015
# Simple test script for Nifti I/O tools

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="nifti_$QUITVER"
mkdir -p $DATADIR
cd $DATADIR

run_test "Create nifti" $QUITDIR/niicreate $DATADIR/blank.nii -d "16 16 16" -v "2 2 2"
run_test "Read header" $QUITDIR/niihdr $DATADIR/blank.nii
SILENCE_TEST="0"
