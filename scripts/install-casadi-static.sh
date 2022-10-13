#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

prefix="$1" # Installation directory prefix
build_type="${2:-RelWithDebInfo}"
version="3.5.5"

if [ -z "$prefix" ]; then
    if [ -z "$VIRTUAL_ENV" ]; then
        echo "No active virtual environment, refusing to install."
        exit 1
    else
        prefix="$VIRTUAL_ENV"
    fi
fi

set -ex
export CMAKE_PREFIX_PATH="$prefix:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$PKG_CONFIG_PATH"

pushd /tmp

# CasADi
[ -d casadi ] \
 || git clone --single-branch --depth=1 --branch "$version" --recursive \
    https://github.com/casadi/casadi
pushd casadi
cmake -S. -Bbuild \
    -G "Ninja Multi-Config" \
    -D CMAKE_INSTALL_PREFIX="$prefix" \
    -D CMAKE_POSITION_INDEPENDENT_CODE=On \
    -D WITH_COMMON=Off \
    -D WITH_PYTHON=Off \
    -D WITH_PYTHON3=Off \
    -D WITH_OPENMP=Off \
    -D WITH_THREAD=On \
    -D WITH_DL=On \
    -D WITH_IPOPT=Off \
    -D ENABLE_STATIC=On \
    -D ENABLE_SHARED=Off \
    -D CMAKE_OSX_ARCHITECTURES="arm64;x86_64"
cmake --build build -j --config $build_type
cmake --install build --config $build_type
popd

popd
