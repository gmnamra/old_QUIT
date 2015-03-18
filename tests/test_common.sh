#!/bin/bash

# Tobias Wood 2015
# Common functions for QUIT Tests

# Simple test function
function run_test {
    "$@"
    local STATUS=$?
    if [ $STATUS -ne 0 ]; then
        echo "Test $1 failed." >&1
        exit $STATUS
    else
		echo "Test $1 passed." >&1
    fi
    return $STATUS
}

# Setup environment
QUITDIR=$PWD/../build
export QUIT_EXT=NIFTI
