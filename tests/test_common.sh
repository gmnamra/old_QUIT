#!/bin/bash

# Tobias Wood 2015
# Common functions for QUIT Tests

# Simple test function
function run_test {
	# $1 is test name, remainder is command to run
	NAME=$1
	shift
	"$@"
	local STATUS=$?
	if [ $STATUS -ne 0 ]; then
		echo "Test $NAME failed." >&1
		exit $STATUS
	else
		echo "Test $NAME passed." >&1
	fi
	return $STATUS
}

# Setup environment
QUITDIR=$PWD/../build
export QUIT_EXT=NIFTI
