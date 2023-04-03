#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..

set -ex

# Determine architecture
case $(uname -m) in
    x86_64) triple="x86_64-centos7-linux-gnu" ;;
    armv6l) triple="armv6-rpi-linux-gnueabihf" ;;
    armv7l) triple="armv7-neon-linux-gnueabihf" ;;
    aarch64) triple="aarch64-rpi3-linux-gnu" ;;
esac

# Download compiler
download_url="https://github.com/tttapa/cross-python/releases/download/0.0.12"
tools_dir="$PWD/toolchains"
pfx="$tools_dir/$triple"
mkdir -p "$tools_dir"
if [ ! -d "$pfx" ]; then
    wget "$download_url/full-$triple.tar.xz" -O- | \
        tar xJ -C "$tools_dir"
fi

# Use ccache to cache compilation
export CMAKE_C_COMPILER_LAUNCHER=ccache
export CMAKE_CXX_COMPILER_LAUNCHER=ccache
# Compile for this system's processor for optimal performance
export CFLAGS="-march=native -fdiagnostics-color"
export CXXFLAGS="-march=native -fdiagnostics-color"
export FCFLAGS="-march=native -fdiagnostics-color"

# Configure
cmake -S. -Bbuild-local \
    --toolchain "$pfx/cmake/$triple.toolchain.cmake" \
    -DCMAKE_FIND_ROOT_PATH="$pfx/eigen;$pfx/casadi;$pfx/openblas;$pfx/mumps;$pfx/ipopt" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=On \
    -DBUILD_SHARED_LIBS=On \
    -DALPAQA_WITH_DRIVERS=On \
    -DALPAQA_WITH_EXAMPLES=Off \
    -DALPAQA_WITH_TESTS=Off \
    -DALPAQA_WITH_PYTHON=Off \
    -DCMAKE_STAGING_PREFIX="$PWD/staging" \
    -G "Ninja Multi-Config"
# Build
for cfg in Debug RelWithDebInfo; do
    cmake --build build-local -j --config $cfg
    cmake --install build-local --config $cfg
    cmake --install build-local --config $cfg --component debug
done
# Package
pushd build-local
cpack -G 'TGZ;DEB' -C "RelWithDebInfo;Debug"
popd

# Build Python package
config="$triple.py-build-cmake.config.toml"
cat <<- EOF > "$config"
[cmake]
config = ["Debug", "Release"]
generator = "Ninja Multi-Config"
[cmake.options]
CMAKE_FIND_ROOT_PATH = "$pfx/pybind11;$pfx/casadi;$pfx/eigen-master"
USE_GLOBAL_PYBIND11 = "On"
ALPAQA_PYTHON_DEBUG_CONFIG = "Debug"
ALPAQA_WITH_PY_STUBS = "On"
EOF
. ./py-venv/bin/activate
LDFLAGS='-static-libgcc -static-libstdc++' \
python -m build -w "." -o staging \
    -C--cross="$pfx/cmake/$triple.py-build-cmake.cross.toml" \
    -C--local="$PWD/$config"
LDFLAGS='-static-libgcc -static-libstdc++' \
python -m build -w "python/alpaqa-debug" -o staging \
    -C--cross="$pfx/cmake/$triple.py-build-cmake.cross.toml" \
    -C--local="$PWD/$config"
pip install -f staging --force-reinstall --no-deps \
    "alpaqa==1.0.0a7" "alpaqa-debug==1.0.0a7"
pip install -f staging \
    "alpaqa[test]==1.0.0a7" "alpaqa-debug==1.0.0a7"
pytest
