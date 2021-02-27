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

# OpenBLAS and LAPACK
[ -d OpenBLAS ] || git clone https://github.com/xianyi/OpenBLAS.git --branch "v0.3.13" --depth 1 --recursive
rm -rf OpenBLAS/build
mkdir -p OpenBLAS/build
pushd OpenBLAS/build
cmake .. \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -D CMAKE_BUILD_TYPE=Release
make -j$(nproc)
make install
popd

popd
