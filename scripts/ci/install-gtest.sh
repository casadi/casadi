#!/usr/bin/env bash

# Script to download and install GoogleTest

set -ex

version="main"
prefix="${1:-$HOME/.local}"

[ -e "$prefix/lib/cmake/GTest/GTestConfig.cmake" ] \
 && exit 0

cd /tmp
# Download
[ -d "googletest" ] \
 || git clone --single-branch --depth=1 --branch "$version" \
        https://github.com/google/googletest.git
mkdir -p googletest/build && cd $_
# Configure
cmake .. -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_BUILD_TYPE=Release
# Build
make -j$(nproc)
# Install
make install
