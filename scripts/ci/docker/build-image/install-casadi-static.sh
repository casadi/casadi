#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

prefix="${1:-/usr/local}"

set -ex
export CMAKE_PREFIX_PATH="$prefix:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$PKG_CONFIG_PATH"

pushd /tmp

# Casadi
[ -d casadi ] || \
    git clone https://github.com/casadi/casadi \
        --branch "3.5.5" --depth 1 --recursive
pushd casadi
cmake -S. -Bbuild \
    -G "Ninja Multi-Config" \
    -D WITH_COMMON=Off \
    -D WITH_PYTHON=Off \
    -D WITH_PYTHON3=Off \
    -D WITH_OPENMP=Off \
    -D WITH_THREAD=On \
    -D WITH_DL=On \
    -D WITH_IPOPT=Off \
    -D CMAKE_INSTALL_PREFIX="$prefix" \
    -D ENABLE_STATIC=On \
    -D ENABLE_SHARED=Off \
    -D CMAKE_POSITION_INDEPENDENT_CODE=On
cmake --build build -j --config Release
cmake --install build --config Release
popd

popd
