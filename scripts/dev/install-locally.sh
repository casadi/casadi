#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..

set -ex

# Download compiler
download_url="https://github.com/tttapa/cross-python/releases/download/0.0.10"
tools_dir="$PWD/toolchains"
triple="x86_64-centos7-linux-gnu"
mkdir -p "$tools_dir"
if [ ! -d "$tools_dir/$triple" ]; then
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
    --toolchain "$tools_dir/$triple/cmake/$triple.toolchain.cmake" \
    -DEigen3_DIR="$tools_dir/$triple/eigen-master/usr/local/share/eigen3/cmake" \
    -Dcasadi_DIR="$tools_dir/$triple/casadi/usr/local/lib/cmake/casadi" \
    -DGTest_DIR="$tools_dir/$triple/googletest/usr/local/lib/cmake/GTest" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=On \
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
staging="$tools_dir/$triple"
config="$triple.py-build-cmake.config.toml"
cat <<- EOF > "$config"
[cmake]
config = ["Debug", "Release"]
generator = "Ninja Multi-Config"
[cmake.options]
CMAKE_FIND_ROOT_PATH = "$staging/pybind11;$staging/casadi;$staging/eigen-master"
USE_GLOBAL_PYBIND11 = "On"
EOF
. ./py-venv/bin/activate
LDFLAGS='-static-libgcc -static-libstdc++' \
python -m build -w "." -o staging \
    -C--cross="$staging/cmake/$triple.py-build-cmake.cross.toml" \
    -C--local="$PWD/$config"
LDFLAGS='-static-libgcc -static-libstdc++' \
python -m build -w "python/alpaqa-debug" -o staging \
    -C--cross="$staging/cmake/$triple.py-build-cmake.cross.toml" \
    -C--local="$PWD/$config"
pip install -f staging --force-reinstall --no-deps \
    "alpaqa==1.0.0a7" "alpaqa-debug==1.0.0a7"
pip install -f staging \
    "alpaqa[test]==1.0.0a7" "alpaqa-debug==1.0.0a7"
pytest
