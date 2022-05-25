#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

prefix="${1:-/usr/local}"
version="3.4.0"

set -ex
export CMAKE_PREFIX_PATH="$prefix:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$PKG_CONFIG_PATH"

pushd /tmp

# Eigen
[ -d eigen ] \
 || git clone --single-branch --depth=1 --branch "$version" \
    https://gitlab.com/libeigen/eigen.git
pushd eigen
cmake -S. -Bbuild \
    -G "Ninja Multi-Config" \
    -D CMAKE_INSTALL_PREFIX="$prefix"
cmake --build build -j --config RelWithDebInfo
cmake --install build --config RelWithDebInfo
popd

popd