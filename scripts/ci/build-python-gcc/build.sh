#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"
set -ex

gcc_version=11
versions=( 3.7 3.8 3.9 3.10 )
for version in "${versions[@]}"; do
    docker build --build-arg PYTHON_VERSION=$version --build-arg GCC_VERSION=$gcc_version -t tttapa/alpaqa-build-python-gcc:$version-$gcc_version .
done
for version in "${versions[@]}"; do
    docker push tttapa/alpaqa-build-python-gcc:$version-$gcc_version
done