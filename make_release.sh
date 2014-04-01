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

set -e
set -x

git tag $1
./update_version.sh
PREFIX=QUIT
ARCHIVE=../${PREFIX}-$1.tar
git archive -v --format tar --prefix ${PREFIX}/ --output ${ARCHIVE} $1
tar -s :src:${PREFIX}/src: -v -r -f ${ARCHIVE} src/version
gzip ${ARCHIVE}
