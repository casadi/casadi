#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

set -ex

pushd /tmp
git clone https://github.com/NixOS/patchelf --branch 0.13
pushd patchelf
./bootstrap.sh
./configure --prefix="$VIRTUAL_ENV"
make -j$(nproc)
make check
make install
popd
popd