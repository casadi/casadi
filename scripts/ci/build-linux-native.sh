#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..
set -ex

# Package and output directories
pkg_dir="${1:-.}"
out_dir="${2:-dist}"
install_stubs_dir="$3"

# Create a py-build-cmake config file
pbc_config="$PWD/native-py-build-cmake.local.pbc"
cat << EOF > "$pbc_config"
conan.profile_host=["default", "$PWD/scripts/ci/profiles/linux-conf.profile"]
conan.cmake.options.CMAKE_C_COMPILER_LAUNCHER=sccache
conan.cmake.options.CMAKE_CXX_COMPILER_LAUNCHER=sccache
conan.cmake.args+=["--fresh"]
conan.cmake.build_args+=["--verbose"]
EOF

# Build the Python package
python3 -m pip install -U build
python3 -m build -w "$pkg_dir" -o "$out_dir" -C local="$pbc_config"

# Install the Python stubs
if [ -n "$install_stubs_dir" ]; then
    # Install py-build-cmake and pybind11-stubgen
    python3 -m pip install 'py-build-cmake~=0.5.1.dev0' 'pybind11-stubgen~=2.5.5' 'numpy<3'
    # Determine Conan's build directory
    pbc=(python3 -m py_build_cmake.cli -C "$pkg_dir" --local="$pbc_config")
    build_config="$("${pbc[@]}" build-config-name)"
    build_dir="$pkg_dir/.py-build-cmake_cache/build/$build_config"
    # Activate the Conan build environment (ensures that CMake is in PATH)
    set +x; source "$build_dir/generators/conanbuild.sh"; set -x
    # Re-run CMake to change Python executable (old one is in a temporary venv)
    cmake "$build_dir" \
        -D "Python3_EXECUTABLE=$(which python3)" \
        -D "Python3_HOST_EXECUTABLE=$(which python3)"
    # Avoid expensive copies of large binary modules
    export CMAKE_INSTALL_MODE=SYMLINK_OR_COPY
    # Install the binary modules into the source tree
    cmake --install "$build_dir" --config Release \
        --prefix "$install_stubs_dir" --component python_modules
    # Generate the stubs (using a complete tree including the binary modules)
    cmake --install "$build_dir" --config Release \
        --prefix "$install_stubs_dir" --component python_stubs
    # Then we remove the binary modules again (sdist is source only)
    while IFS= read -r f || [ -n "$f" ]; do rm -f "$f"
    done < "$build_dir/install_manifest_python_modules.txt"
fi
