#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

./scripts/install-openblas.sh
./scripts/install-eigen.sh
./scripts/install-casadi.sh
./scripts/install-cutest.sh
./scripts/install-yaml-cpp.sh
./scripts/install-gtest.sh
