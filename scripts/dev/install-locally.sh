#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..
export CTEST_OUTPUT_ON_FAILURE=1

set -ex

build_python=1
build_cpp=1

case $1 in
    "") ;;
    py*) build_cpp=0 ;;
    cpp*) build_python=0 ;;
    *) echo "Invalid argument, use py[-dev] or cpp[-pkg]"; exit 1 ;;
esac

# Activate virtual environment
[ -d ".venv" ] || python3 -m venv .venv
. ./.venv/bin/activate
python_version="$(python --version | cut -d ' ' -f 2)"
python_majmin="${python_version%.*}"
python_majmin_nodot="${python_majmin//./}"

# Select architecture
triple=x86_64-bionic-linux-gnu
plat_tag=manylinux_2_27_x86_64

# Download dependencies
pip install -U pip build conan
# My own custom recipes for Ipopt, CasADi, QPALM, patched Eigen
tools_dir="$PWD/toolchains"
[ -d "$tools_dir/thirdparty/conan-recipes" ] || {
    mkdir -p "$tools_dir/thirdparty"
    git clone https://github.com/tttapa/conan-recipes "$tools_dir/thirdparty/conan-recipes"
    conan remote add tttapa-conan-recipes "$tools_dir/thirdparty/conan-recipes" --force
}

# Create Conan profiles
host_profile="$PWD/profile-host.local.conan"
cat <<- EOF > "$host_profile"
include($PWD/scripts/ci/profiles/$triple.profile)
[conf]
tools.build.cross_building:can_run=True
[buildenv]
CMAKE_C_COMPILER_LAUNCHER=sccache
CMAKE_CXX_COMPILER_LAUNCHER=sccache
EOF

cpp_profile="$PWD/profile-cpp.local.conan"
cat <<- EOF > "$cpp_profile"
include($host_profile)
include($PWD/scripts/ci/profiles/alpaqa-cpp-linux.profile)
EOF

python_profile="$PWD/profile-python.local.conan"
cat <<- EOF > "$python_profile"
include($host_profile)
include($PWD/scripts/ci/profiles/alpaqa-python-linux.profile)
[conf]
tools.build:exelinkflags+=["-static-libgcc"]
tools.build:sharedlinkflags+=["-static-libgcc"]
[options]
alpaqa/*:with_conan_python=True
[replace_requires]
tttapa-python-dev/*: tttapa-python-dev/[^$python_version]
EOF

pbc_config="$PWD/$triple.py-build-cmake.config.pbc"
cat <<- EOF > "$pbc_config"
os=linux
implementation="cp"
version="$python_majmin_nodot"
abi="cp$python_majmin_nodot"
arch=$plat_tag
cmake.options.CMAKE_C_COMPILER_LAUNCHER=sccache
cmake.options.CMAKE_CXX_COMPILER_LAUNCHER=sccache
cmake.options.ALPAQA_WITH_PY_STUBS=true
EOF

# Build C++ packages
if [ $build_cpp -eq 1 ]; then
    for cfg in Debug RelWithDebInfo; do
        conan install . --build=missing \
            -pr:h "$cpp_profile" \
            -s build_type=$cfg
    done

    # Configure
    cmake --preset conan-default -B build \
        -D CMAKE_EXPORT_COMPILE_COMMANDS=On \
        -D CMAKE_INSTALL_PREFIX="$PWD/staging"
    # Build
    for cfg in Debug RelWithDebInfo; do
        cmake --build build -j --config $cfg
        cmake --install build --config $cfg
        cmake --install build --config $cfg --component debug
    done
    # Package
    if [ "${1: -4}" = "-pkg" ]; then
        pushd build
        cpack -G 'TGZ;DEB' -C "RelWithDebInfo;Debug"
        popd
    fi
fi

# Build Python packages
if [ $build_python -eq 1 ]; then
    for cfg in Debug Release; do
        conan install . --build=missing \
            -pr:h "$python_profile" \
            -s build_type=$cfg
    done
    if [ "${1: -4}" = "-dev" ]; then
        pip install -e ".[test]" -v \
            -C local="$PWD/scripts/ci/py-build-cmake.toml" \
            -C cross="$pbc_config" \
            -C override=cmake.install_components+='[python_modules_debug]'
    else
        tag=\"$(date -u +"%s")\"
        python -m build -w "." -o staging \
            -C local="$PWD/scripts/ci/py-build-cmake.toml" \
            -C cross="$pbc_config" \
            -C override=wheel.build_tag="$tag"
        python -m build -w "python/alpaqa-debug" -o staging \
            -C local="$PWD/scripts/ci/py-build-cmake.toml" \
            -C component="$PWD/scripts/ci/py-build-cmake.component.toml" \
            -C cross="$pbc_config" \
            -C override=wheel.build_tag="$tag"
        pip install -f staging --force-reinstall --no-deps \
            "alpaqa==1.0.0a21.dev0" "alpaqa-debug==1.0.0a21.dev0"
        pip install -f staging \
            "alpaqa[test]==1.0.0a21.dev0" "alpaqa-debug==1.0.0a21.dev0"
    fi
    pytest
fi
