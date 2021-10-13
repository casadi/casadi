#!/usr/bin/env bash

# Script to download and install Eigen

set -ex

version="3.4.0"
prefix="${1:-$HOME/.local}"

[ -e "$prefix/share/eigen3/cmake/Eigen3Config.cmake" ] \
 && exit 0

cd /tmp
# Download
[ -d "eigen" ] \
 || git clone --single-branch --depth=1 --branch "$version" \
    https://gitlab.com/libeigen/eigen.git
mkdir -p eigen/build && cd $_
# Configure
cmake .. -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_BUILD_TYPE=Release
# Install
make install
