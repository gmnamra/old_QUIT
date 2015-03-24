#!/bin/bash -e

# make_release.sh
# DESPOT
#
# Created by Tobias Wood on 05/03/2014.
#
# Merges development branch to master, sets the version string
# and then applies the same string as a tag, then pushes the
# whole lot to github.

TAG=$(date +%y%m%d)
VERSION="\"QUIT_$TAG\""

STATUS=$(git status -s -uno)

if [ -n "$STATUS" ]; then
	echo "YOU HAVE UNCOMMITED CHANGES. ABORTING."
	exit 1
fi

git checkout master
# Use no-ff to make it clear we are merging a branch
git merge --no-ff development -m "Merging development branch for release $TAG"
echo $VERSION > src/version
git commit src/version -m "Updated version file for release $TAG"
git tag $TAG
git push github master --tags

git checkout development
