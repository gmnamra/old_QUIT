#!/bin/bash

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh

DATADIR="despot"
mkdir -p $DATADIR
cd $DATADIR

echo "Starting DESPOT tests"

DIMS="32 32 32"
run_test "Create T1 volume" $QUITDIR/niicreate -d "$DIMS" -g "0 0.5 5" T1.nii
run_test "Create T2 volume" $QUITDIR/niicreate -d "$DIMS" -g "1 0.05 0.5" T2.nii
run_test "Create PD volume" $QUITDIR/niicreate -d "$DIMS" -f 1 PD.nii
run_test "Create f0 volume" $QUITDIR/niicreate -d "$DIMS" -g "2 -50.0 50.0" f0.nii

# Create input for mcsignal
INPUT="PD.nii
T1.nii
T2.nii
f0.nii
SPGR
spgr.nii
3
5 10 15
0.01
SSFP
ssfp.nii
3
15 30 45
4
0 90 180 270
0.005
END"
echo "$INPUT" > mcsignal.in
run_test "Create SPGR data" $QUITDIR/mcsignal --1 < mcsignal.in
cd ..
