#!/bin/bash

#  release.sh
#  DESPOT
#
#  Created by Tobias Wood on 05/03/2014.
#

if [ -z "$1" ]
then
	echo "Must supply a release name."
	exit 1
fi

git tag $1
./update_version.sh
PREFIX=QUIT
git archive -v --format tar.gz --prefix ${PREFIX}/ --output ../${PREFIX}-$1.tar.gz $1
