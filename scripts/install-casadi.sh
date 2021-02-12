#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

set -ex
python --version
which python

pushd /tmp

# SWIG
rm -rf swig
git clone https://github.com/swig/swig --branch v4.0.2 --depth 1
pushd swig 
./autogen.sh 
./configure --prefix="$VIRTUAL_ENV"
make -j$(nproc)
make install
export CMAKE_PREFIX_PATH="$VIRTUAL_ENV:$CMAKE_PREFIX_PATH"
popd

# Casadi
[ -d casadi ] || git clone https://github.com/casadi/casadi --branch 3.5.5 --depth 1 --recursive
rm -rf casadi/build
mkdir "$_"
pushd "$_"
cmake .. \
    -D WITH_COMMON=On \
    -D WITH_PYTHON=On \
    -D WITH_PYTHON3=On \
    -D WITH_OPENMP=On \
    -D WITH_THREAD=On \
    -D WITH_DL=On \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo
make -j$(nproc)
make install
popd

popd
