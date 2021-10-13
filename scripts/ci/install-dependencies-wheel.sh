#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/../..

set -ex

# ./scripts/install-openblas.sh Release
./scripts/install-eigen.sh Release
./scripts/ci/install-casadi-wheel.sh Release
# ./scripts/install-cutest.sh
./scripts/install-yaml-cpp.sh Release
./scripts/install-gtest.sh
./scripts/install-lbfgspp.sh Release