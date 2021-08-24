#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"
set -ex

versions=( 3.7 3.8 3.9 rc )
for version in "${versions[@]}"; do
    docker build --build-arg PYTHON_VERSION=$version -t tttapa/build-python-gcc:$version-10 .
done
for version in "${versions[@]}"; do
    docker push tttapa/build-python-gcc:$version-10
done