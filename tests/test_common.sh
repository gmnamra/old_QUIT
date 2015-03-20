#!/bin/bash

# Tobias Wood 2015
# Common functions for QUIT Tests

# Check for presence of FSL
if [ -n "$(which fslmaths)" ]; then
	HAVE_FSL="1"
	echo "FSL detected."
else
	HAVE_FSL=""
	echo "FSL not detected. Cannot run comparison tests."
fi

# Simple test function
SILENCE_TESTS=0
function run_test {
	# $1 is test name, remainder is command to run
	NAME="$1"
	shift
	echo "Starting test $NAME"
	if [ "$SILENCE_TESTS" -eq "1" ]; then
		"$@" > "$NAME.log"
	else
		"$@"
	fi
	local STATUS=$?
	if [ $STATUS -ne 0 ]; then
		echo "Test $NAME failed." >&2
		exit $STATUS
	else
		echo "Test $NAME passed." >&1
	fi
	return $STATUS
}

# Use FSL to compare 2 images and check if the average difference is below the tolerance
function compare_test {
	# $1 test name
	# $2 reference, $3 test image, $4 tolerance
	NAME="$1"
	REF=$2
	TEST=$3
	TOL=$4
	DIFF=${REF%.nii}_${TEST%.nii}
	if [ "$HAVE_FSL" -eq "1" ]; then
		fslmaths $REF -sub $TEST $DIFF
		MEAN=$(fslstats $DIFF -M)
		TEST=$(echo "$MEAN $TOL" | awk ' { if(sqrt($1*$1)<=$2) { print 1 } else { print 0 }}')
		if [ "$TEST" -eq "1" ]; then
			echo "Comparison test $NAME passed, value was $MEAN tolerance $TOL"
		else
			echo "Comparison test $NAME failed, value was $MEAN tolerance $TOL"
			exit 1
		fi
	else
		echo "FSL not present, skipping test $NAME"
	fi
	return 0
}

# Setup environment
QUITDIR=$PWD/../build
export QUIT_EXT=NIFTI
QUITVER=$(cat ../src/version | sed -e 's/^"//'  -e 's/"$//')
