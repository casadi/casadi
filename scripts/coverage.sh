#!/usr/bin/env bash

set -ex

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

proj_dir=$(dirname "$dir")
build_dir=$(pwd)

html_dest="$proj_dir/docs/Coverage"
dest="$build_dir/coverage"
mkdir -p "$dest"
mkdir -p "$html_dest"

rm -f "$dest/*.info"
rm -rf "$html_dest/*"

# Parse command line arguments

compiler="${1,,}"
version="${2%%.*}"
compiler_path="${3}"

if [[ "$compiler_path" == *"-g++" ]]; then
    GCOV_BIN="${compiler_path%-g\+\+}-gcov"
    CPPFILT_BIN="${compiler_path%-g\+\+}-c++filt"
fi

if [ -n "$version" ]; then version="-${version}"; fi

echo "Compiler: ${compiler}${version} (${compiler_path})"

# If the compiler is Clang, use a wrapper around llvm-cov that emulates gcov
if [ -n "$GCOV_BIN" ]; then
    gcov_tool=("--gcov-tool" "$GCOV_BIN");
elif [ "${compiler}" == "clang" ]; then
    gcov_tool=("--gcov-tool" "llvm-cov${version}" "--gcov-tool" "gcov")
else
    gcov_tool=("--gcov-tool" "gcov${version}")
fi

# Replace the default c++filt program with LLVM/Clang's version
if [ -n "$CPPFILT_BIN" ]; then
    cppfilt_bin="$CPPFILT_BIN"
elif [ "${compiler}" == "clang" ]; then
    cppfilt_bin="$(which llvm-cxxfilt${version})"
else
    cppfilt_bin="$(which c++filt)"
fi

branches=0

# Reset counters
lcov \
    --zerocounters \
    --directory "$build_dir"

# Initial capture
lcov \
    --capture --initial \
    --directory "$build_dir" \
    --include "$proj_dir"'/src/alpaqa/**' \
    --include "$proj_dir"'/src/interop/**' \
    --output-file "$dest"/coverage_base.info \
    "${gcov_tool[@]}" \
    --rc branch_coverage=$branches

# Run tests
ctest

# Actual capture
lcov \
    --capture \
    --directory "$build_dir" \
    --include "$proj_dir"'/src/alpaqa/**' \
    --include "$proj_dir"'/src/interop/**' \
    --output-file "$dest"/coverage_test.info \
    "${gcov_tool[@]}" \
    --rc branch_coverage=$branches \
    --rc check_data_consistency=0

# Combine captures
lcov \
    --add-tracefile "$dest"/coverage_base.info \
    --add-tracefile "$dest"/coverage_test.info \
    --output-file "$dest"/coverage_total.info \
    "${gcov_tool[@]}" \
    --rc branch_coverage=$branches \
    --ignore-errors corrupt,inconsistent

# Generate HTML coverage report
genhtml \
    --prefix "$proj_dir" \
    "$dest"/coverage_total.info \
    --output-directory="$html_dest" \
    --legend --title $(cd "$proj_dir" && git rev-parse HEAD) \
    --rc branch_coverage=$branches \
    -s --demangle-cpp "$cppfilt_bin" \
    --no-function-coverage # because of the many templates

# python3 "$proj_dir/scripts/coverage-badge.py"
