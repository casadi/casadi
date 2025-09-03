#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..
export CTEST_OUTPUT_ON_FAILURE=1

set -ex

# Select architecture
triple="${1:-x86_64-bionic-linux-gnu}"
tests="${2}"

if [ -n "$tests" ]; then test_flags="-DALPAQA_FORCE_TEST_DISCOVERY=On"; fi

# Create Conan profiles
cpp_profile="$PWD/profile-cpp.conan"
cat <<- EOF > "$cpp_profile"
include($PWD/scripts/ci/profiles/$triple.profile)
include($PWD/scripts/ci/profiles/alpaqa-cpp-linux.profile)
[conf]
tools.cmake.cmaketoolchain:user_toolchain=+['$PWD/scripts/ci/profiles/static-libgcc.cmake']
tools.cmake.cmake_layout:build_folder_vars=['const.pkg']
*:tools.build:skip_test=True
&:tools.build:skip_test=False
[buildenv]
LDFLAGS+= -static-libstdc++ -static-libgfortran -static-libquadmath -Wl,--as-needed 
&:CMAKE_C_COMPILER_LAUNCHER=sccache
&:CMAKE_CXX_COMPILER_LAUNCHER=sccache
[options]
coinmumps/*:static_fortran_libs=True
[replace_requires]
eigen/*: eigen/3.4.0
EOF

for cfg in Debug Release; do
    # Dependencies
    conan install . --build=missing -pr:h "$cpp_profile" -s build_type=$cfg
    . ./build/pkg/generators/conanbuild.sh
    # Configure
    cmake --preset conan-pkg --fresh $test_flags
    # Build
    cmake --build build/pkg -j --config $cfg
    # Test
    if [ -n "$tests" ]; then
        ctest -C $cfg --test-dir build/pkg -D $tests
    fi
    . ./build/pkg/generators/deactivate_conanbuild.sh
done
# Package
pushd build/pkg
. ./generators/conanbuild.sh
cpack -G 'TGZ;DEB' -C "Release;Debug"
popd
