#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..

set -ex

# Matlab folder
matlab_dir="${1:-/usr/local/MATLAB}"

# Select architecture
triple="${2:-x86_64-bionic-linux-gnu}"

# Package and output directories
pkg_dir="${3:-.}"
out_dir="${4:-staging/matlab}"

# Create Conan profiles
matlab_profile="$PWD/profile-matlab.local.conan"
cat <<- EOF > "$matlab_profile"
include($PWD/scripts/ci/profiles/$triple.profile)
include($PWD/scripts/ci/profiles/alpaqa-matlab.profile)
[conf]
tools.build:exelinkflags+=["-static-libgcc"]
tools.build:sharedlinkflags+=["-static-libgcc"]
EOF

# Install dependencies
rm -rf "$pkg_dir"/build/matlab-{debug,release}/{generators,CMakeCache.txt}
for cfg in Release; do
    conan install "$pkg_dir" --build=missing \
        -pr:h "$matlab_profile" \
        -s build_type=$cfg
done

# Build MATLAB bindings
pushd "$pkg_dir"
cmake --preset conan-matlab-release \
    -D CMAKE_FIND_ROOT_PATH="$matlab_dir" \
    -D CMAKE_C_COMPILER_LAUNCHER=sccache \
    -D CMAKE_CXX_COMPILER_LAUNCHER=sccache
cmake --build --preset conan-matlab-release \
    -t alpaqa_mex -v
cmake --install build/matlab-release \
    --component mex_interface --prefix "$out_dir"
popd
