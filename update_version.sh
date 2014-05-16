#!/bin/bash

#  update_version.sh
#  DESPOT
#
#  Created by Tobias Wood on 05/03/2014.
#

VER="const std::string version { \"QUIT `git describe --tags`\" };"
echo "$VER" > src/version
