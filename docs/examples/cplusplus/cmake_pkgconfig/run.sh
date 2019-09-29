#!/bin/bash

#export PKG_CONFIG_PATH=..
rm -rf build
mkdir build
pushd build
cmake ..
make VERBOSE=1
popd
build/casadi_demo


