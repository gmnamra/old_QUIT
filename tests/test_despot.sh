#!/bin/bash

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

DATADIR="data"
QUITDIR=$PWD/../build
mkdir -p $DATADIR
cd $DATADIR

$QUITDIR/niicreate

cd ..
