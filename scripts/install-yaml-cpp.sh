#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

set -ex
export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"

pushd /tmp

# yaml-cpp
[ -d yaml-cpp ] || git clone https://github.com/jbeder/yaml-cpp.git --branch yaml-cpp-0.6.3 --depth 1 --recursive
rm -rf yaml-cpp/build
mkdir -p yaml-cpp/build
pushd yaml-cpp/build
cmake .. \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -D YAML_BUILD_SHARED_LIBS=On \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo \
    -D YAML_CPP_INSTALL=On
make -j$(nproc)
make install
popd

popd
