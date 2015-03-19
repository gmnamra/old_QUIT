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

DIMS="16 16 16"
run_test "Create T1 volume" $QUITDIR/niicreate -d "$DIMS" -f 1 T1.nii
run_test "Create T2 volume" $QUITDIR/niicreate -d "$DIMS" -g "2 0.05 0.5" T2.nii
run_test "Create PD volume" $QUITDIR/niicreate -d "$DIMS" -f 1 PD.nii

run_test "Create SPGR data" $QUITDIR/mcsignal
cd ..
