#!/bin/bash

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="mcdespot_$QUITVER"
mkdir -p $DATADIR
cd $DATADIR

DIMS="11 21 1"
#DIMS="3 3 3"

$QUITDIR/niicreate -d "$DIMS" -f "1.0" PD.nii
$QUITDIR/niicreate -d "$DIMS" -f "0.25" T1_a.nii
$QUITDIR/niicreate -d "$DIMS" -f "0.015" T2_a.nii
$QUITDIR/niicreate -d "$DIMS" -f "1.5" T1_b.nii
$QUITDIR/niicreate -d "$DIMS" -f "0.1" T2_b.nii
$QUITDIR/niicreate -d "$DIMS" -g "0 0 1.0" f_a.nii
$QUITDIR/niicreate -d "$DIMS" -f "0.1" tau_a.nii
$QUITDIR/niicreate -d "$DIMS" -g "1 -25. 25." f0.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_PAR="8
3 5 8 11 14 17 20 23
0.01"
SSFP_FILE="ssfp.nii"
SSFP_PAR="8
5 10 20 30 40 50 60 70
2
180 0
0.05"

# Create input for Multi-Component
MCSIG_INPUT="PD.nii
T1_a.nii
T2_a.nii
T1_b.nii
T2_b.nii
tau_a.nii
f_a.nii
f0.nii
SPGR
$SPGR_FILE
$SPGR_PAR
SSFP
$SSFP_FILE
$SSFP_PAR
END"
echo "$MCSIG_INPUT" > mcsignal.in
run_test "CREATE_SIGNALS" $QUITDIR/mcsignal --2 < mcsignal.in

echo "SPGR
$SPGR_FILE
$SPGR_PAR
SSFP
$SSFP_FILE
$SSFP_PAR" > mcdespot.in

run_test "MCDESPOT" $QUITDIR/mcdespot --2 -n -v -S 1 -T 1 < mcdespot.in
compare_test "MWF" f_a.nii 2C_f_a.nii 0.1

cd ..
SILENCE_TESTS="0"
