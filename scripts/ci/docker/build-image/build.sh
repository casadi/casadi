#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"
set -ex

gcc_version=11
docker buildx build --build-arg GCC_VERSION=$gcc_version \
    --target alpaqa-build -t tttapa/alpaqa-build .
docker push tttapa/alpaqa-build &
docker buildx build --build-arg GCC_VERSION=$gcc_version \
    --target alpaqa-test -t tttapa/alpaqa-test .
docker push tttapa/alpaqa-test &
versions=( 3.8 3.9 3.10 )
for version in "${versions[@]}"; do
    docker buildx build --build-arg GCC_VERSION=$gcc_version --build-arg PYTHON_VERSION=$version \
        --target alpaqa-python-build -t tttapa/alpaqa-python-build:py$version .
    docker push tttapa/alpaqa-python-build:py$version &
done
wait