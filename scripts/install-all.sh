#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

rm -rf py-venv
python3.9 -m venv py-venv
. ./py-venv/bin/activate
pip install -r scripts/requirements.txt

./scripts/install-openblas.sh   # https://www.openblas.net/
./scripts/install-eigen.sh      # https://eigen.tuxfamily.org/index.php
./scripts/install-casadi.sh     # https://web.casadi.org/
./scripts/install-cutest.sh     # https://github.com/ralna/CUTEst
./scripts/install-lbfgspp.sh    # https://github.com/yixuan/LBFGSpp
./scripts/install-yaml-cpp.sh   # https://github.com/jbeder/yaml-cpp
./scripts/install-gtest.sh      # https://google.github.io/googletest/

mkdir build
pushd build
cmake .. -D CMAKE_BUILD_TYPE=Release
make -j$(nproc) casadi-rosenbrock cutest-rosenbrock
./examples/CasADi/Rosenbrock/casadi-rosenbrock
./examples/CUTEst/Rosenbrock/cutest-rosenbrock
popd
