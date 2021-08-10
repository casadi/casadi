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

# LBFGS++
[ -d LBFGSpp ] || \
git clone --single-branch --depth=1 --branch master \
    https://github.com/tttapa/LBFGSpp.git
rm -rf LBFGSpp/build
mkdir -p LBFGSpp/build
pushd LBFGSpp/build
cmake .. \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -D CMAKE_BUILD_TYPE="${build_type}"
make -j$(nproc)
make install
popd

popd
