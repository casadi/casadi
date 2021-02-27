#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

rm -rf py-venv
python3.9 -m venv py-venv
. ./py-venv/bin/activate
pip install -r scripts/requirements.txt

./scripts/install-openblas.sh
./scripts/install-eigen.sh
./scripts/install-casadi.sh
./scripts/install-cutest.sh
./scripts/install-yaml-cpp.sh
./scripts/install-gtest.sh

mkdir build
pushd build
cmake .. -D CMAKE_BUILD_TYPE=Release
make -j8 casadi-rosenbrock cutest-rosenbrock
./examples/CasADi/Rosenbrock/casadi-rosenbrock
./examples/CUTEst/Rosenbrock/cutest-rosenbrock
popd
