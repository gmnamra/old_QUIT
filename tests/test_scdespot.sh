#!/bin/bash

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="scdespot_$QUITVER"
mkdir -p $DATADIR
cd $DATADIR

DIMS="11 11 11"
$QUITDIR/niicreate -d "$DIMS" -f 1 PD.nii
$QUITDIR/niicreate -d "$DIMS" -g "0 0.5 5" T1.nii
$QUITDIR/niicreate -d "$DIMS" -g "1 0.05 0.5" T2.nii
$QUITDIR/niicreate -d "$DIMS" -g "2 -25.0 25.0" f0.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_PAR="3
5 10 15
0.01"
SSFP_FILE="ssfp.nii"
SSFP_PAR="3
15 30 45
4
0 90 180 270
0.05"
MPRAGE_FILE="mprage.nii"
MPRAGE_PAR="5
0.01
11
1
0.5
0"

# Create input for Single Component
MCSIG_INPUT="PD.nii
T1.nii
T2.nii
f0.nii
SPGR
$SPGR_PAR
$SPGR_FILE
SSFP
$SSFP_PAR
$SSFP_FILE
MPRAGE
$MPRAGE_PAR
$MPRAGE_FILE
END"
echo "$MCSIG_INPUT" > mcsignal.in
run_test "CREATE_SIGNALS" $QUITDIR/mcsignal --1 < mcsignal.in

echo "$SPGR_PAR" > despot1.in
echo "$SPGR_PAR
$MPRAGE_PAR" > despot1hifi.in
echo "$SSFP_PAR" > despot2fm.in
D2GS_PAR="3
15 30 45
0.05"
echo "$D2GS_PAR" > despot2gs.in

run_test "DESPOT1" $QUITDIR/despot1 $SPGR_FILE -n < despot1.in
run_test "DESPOT1LM" $QUITDIR/despot1 $SPGR_FILE -n -an -oN < despot1.in
run_test "DESPOT1HIFI" $QUITDIR/despot1hifi $SPGR_FILE $MPRAGE_FILE -n -N < despot1hifi.in
run_test "SSFPBANDS" $QUITDIR/ssfpbands -22 $SSFP_FILE
run_test "DESPOT2FM" $QUITDIR/despot2fm D1_T1.nii $SSFP_FILE -n -fASYM -v < despot2fm.in
run_test "DESPOT2GS" $QUITDIR/despot2 -e D1_T1.nii ${SSFP_FILE%.nii}_lreg_2p.nii -n < despot2gs.in
compare_test "DESPOT1" T1.nii D1_T1.nii 0.01
compare_test "DESPOT1LM" T1.nii ND1_T1.nii 0.01
compare_test "HIFI_T1" T1.nii HIFI_T1.nii 0.01
compare_test "DESPOT2FM" T2.nii FM_T2.nii 0.05
compare_test "DESPOT2GS" T2.nii D2_T2.nii 0.05

cd ..
SILENCE_TESTS="0"
