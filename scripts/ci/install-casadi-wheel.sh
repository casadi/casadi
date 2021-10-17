#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

build_type="${1:-RelWithDebInfo}"

set -ex
export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$VIRTUAL_ENV/lib/pkgconfig:$PKG_CONFIG_PATH"

pushd /tmp

# Casadi
[ -d casadi ] || git clone https://github.com/casadi/casadi --branch "3.5.5" --depth 1 --recursive
rm -rf casadi/build
pushd casadi
cmake -Bbuild -S. -D CMAKE_BUILD_TYPE="${build_type}" \
    -D WITH_COMMON=Off \
    -D WITH_OPENMP=Off \
    -D WITH_THREAD=On \
    -D WITH_DL=On \
    -D ENABLE_STATIC=On \
    -D ENABLE_SHARED=Off
cmake --build build -j1
cmake --install build
popd

popd
