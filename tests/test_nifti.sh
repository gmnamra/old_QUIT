#!/bin/bash
set -ex

# Tobias Wood 2015
# Simple test script for Nifti I/O tools

# Create some test files

QUITDIR=$PWD/../build
export QUIT_EXT=NIFTI
DATADIR="nii_data"
mkdir -p $DATADIR

$QUITDIR/niicreate $DATADIR/blank.nii 16 16 16 1 1 1 1 1
$QUITDIR/niihdr $DATADIR/blank.nii

echo "Tests complete"
