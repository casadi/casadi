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

# SWIG
rm -rf swig
git clone https://github.com/swig/swig --branch "v4.0.2" --depth 1
pushd swig
./autogen.sh
./configure --prefix="$VIRTUAL_ENV"
make -j$(nproc)
make install
popd

# Casadi
[ -d casadi ] || git clone https://github.com/casadi/casadi --branch "3.5.5" --depth 1 --recursive
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
    -D WITH_IPOPT=Off \
    -D WITH_MUMPS=Off \
    -D ENABLE_STATIC=On \
    -D ENABLE_SHARED=Off \
    -D CMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -D CMAKE_BUILD_TYPE="${build_type}" \
    -D PYTHON_VERSION_STRING="$(python --version | cut -f2 -d' ')" \
    -D PYTHON_EXECUTABLE="$(which python)" \
    -D PYTHON_LIBRARY="$(python -c 'import sysconfig; print(sysconfig.get_config_var("LIBDIR") + "/" + sysconfig.get_config_var("LDLIBRARY"))')" \
    -D PYTHON_INCLUDE_DIR="$(python -c 'import sysconfig; print(sysconfig.get_path("include"))')"
make -j$(nproc)
make install
popd

popd
