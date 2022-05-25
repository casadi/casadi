#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

prefix="${1:-/usr/local}"

set -ex
export CMAKE_PREFIX_PATH="$prefix:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$PKG_CONFIG_PATH"

pushd /tmp

# Google Test
[ -d googletest ] \
 || git clone --single-branch --depth=1 --branch main \
    https://github.com/google/googletest.git
pushd googletest
cmake -S. -Bbuild \
    -G "Ninja Multi-Config" \
    -D CMAKE_INSTALL_PREFIX="$prefix"
cmake --build build -j --config Release
cmake --install build --config Release
cmake --build build -j --config Debug
cmake --install build --config Debug
popd

popd