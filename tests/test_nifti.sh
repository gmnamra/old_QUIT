#!/bin/bash

# Tobias Wood 2015
# Simple test script for Nifti I/O tools

source ./test_common.sh

DATADIR="nii_data"
mkdir -p $DATADIR

echo "Starting NIFTI tests."

run_test "Create nifti" $QUITDIR/niicreate $DATADIR/blank.nii -d "16 16 16" -v "2 2 2"
run_test "Read header" $QUITDIR/niihdr $DATADIR/blank.nii
