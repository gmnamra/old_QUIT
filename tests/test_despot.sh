#!/bin/bash

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh

DATADIR="data"
mkdir -p $DATADIR
cd $DATADIR

echo "Starting DESPOT tests"

run_test $QUITDIR/niicreate

cd ..
