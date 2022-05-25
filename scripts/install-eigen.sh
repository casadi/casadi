#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

build_type="${1:-RelWithDebInfo}"
version="3.4.0"

set -ex
export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$VIRTUAL_ENV/lib/pkgconfig:$PKG_CONFIG_PATH"

pushd /tmp

# Eigen
[ -d eigen ] \
 || git clone --single-branch --depth=1 --branch "$version" \
    https://gitlab.com/libeigen/eigen.git
pushd eigen
cmake -S. -Bbuild \
    -G "Ninja Multi-Config" \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV"
cmake --build build -j --config RelWithDebInfo
cmake --install build --config RelWithDebInfo
popd

popd
