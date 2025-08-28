#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..
set -ex

# Package and output directories
pkg_dir="${1:-.}"
out_dir="${2:-dist}"
build_type="${3:-debug}"  # debug, release or relwithdebinfo

# Create a py-build-cmake config file
pbc_config="$PWD/native-py-build-cmake.local.pbc"
cat << EOF > "$pbc_config"
conan.profile_host=["default", "$PWD/scripts/ci/profiles/linux-conf-hardened-$build_type.profile"]
conan.cmake.options.CMAKE_C_COMPILER_LAUNCHER=sccache
conan.cmake.options.CMAKE_CXX_COMPILER_LAUNCHER=sccache
conan.cmake.args+=["--fresh"]
conan.cmake.build_args+=["--verbose"]
conan.cmake.build_type=!
EOF

# Build the Python package
python3 -m pip install -U build
python3 -m build -w "$pkg_dir" -o "$out_dir" -C local="$pbc_config"
