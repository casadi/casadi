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

# MUMPS linear solver
rm -rf ThirdParty-Mumps
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git --branch "master" --depth 1 --recursive
pushd ThirdParty-Mumps
./get.Mumps
./configure --prefix="$VIRTUAL_ENV" --with-lapack-lflags="$(pkg-config --libs-only-l openblas) -pthread -lm" LDFLAGS="$(pkg-config --libs-only-L openblas)"
make -j$(nproc)
make install
popd

# Ipopt
[ -d Ipopt ] || git clone https://github.com/coin-or/Ipopt.git --branch "releases/3.13.3" --depth 1 --recursive
rm -rf Ipopt/build
mkdir -p Ipopt/build
pushd Ipopt/build
../configure --prefix="$VIRTUAL_ENV" --with-lapack-lflags="$(pkg-config --libs-only-l openblas) -pthread -lm" LDFLAGS="$(pkg-config --libs-only-L openblas)"
make -j$(nproc)
make install
popd

# SWIG
rm -rf swig
git clone https://github.com/swig/swig --branch v4.0.2 --depth 1
pushd swig 
./autogen.sh 
./configure --prefix="$VIRTUAL_ENV"
make -j$(nproc)
make install
popd

# Casadi
[ -d casadi ] || git clone https://github.com/casadi/casadi --branch 3.5.5 --depth 1 --recursive
rm -rf casadi/build
mkdir -p casadi/build
pushd casadi/build
cmake .. \
    -D WITH_COMMON=On \
    -D WITH_PYTHON=On \
    -D WITH_PYTHON3=On \
    -D WITH_OPENMP=On \
    -D WITH_THREAD=On \
    -D WITH_DL=On \
    -D WITH_IPOPT=On \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo
make -j$(nproc)
make install
popd

popd
