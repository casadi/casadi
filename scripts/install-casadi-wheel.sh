#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

set -ex
export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$VIRTUAL_ENV/lib/pkgconfig:$PKG_CONFIG_PATH"

pushd /tmp

# CasADi
[ -d casadi ] || \
    git clone https://github.com/casadi/casadi \
        --branch "3.5.5" --depth 1 --recursive
pushd casadi
cmake -Bbuild -S. \
    -G "Ninja Multi-Config" \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -D WITH_COMMON=Off \
    -D WITH_OPENMP=Off \
    -D WITH_THREAD=On \
    -D WITH_DL=On \
    -D ENABLE_STATIC=On \
    -D ENABLE_SHARED=Off
cmake --build build -j --config Release
cmake --install build --config Release
popd

popd
