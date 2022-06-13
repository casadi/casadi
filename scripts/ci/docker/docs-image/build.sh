#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"
set -ex

docker buildx build \
    --target alpaqa-docs -t tttapa/alpaqa-docs .
docker push tttapa/alpaqa-docs &
wait