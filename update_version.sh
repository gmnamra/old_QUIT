#!/bin/bash

#  update_version.sh
#  DESPOT
#
#  Created by Tobias Wood on 05/03/2014.
#

VER="const std::string version { \"DESPOT Tools `git describe --tags`\" };"
echo "$VER" > src/version
