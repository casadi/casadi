# Installation for developers #{pages-installation-dev}

## Standard dependencies

- Python 3.9
- GCC 10 or Clang 10
- Gfortran 10
- CMake
- GNU Make
- GNU Autotools and Automake
- PCRE3

On Ubuntu:
```sh
sudo apt install python3 gcc-10 g++-10 gfortran-10 cmake make autotools-dev automake libpcre3-dev
```
Make sure that a default version of GCC is selected:
```sh
sudo update-alternatives --remove-all gcc ||:
sudo update-alternatives --remove-all g++ ||:
sudo update-alternatives --remove-all gfortran ||:

sudo update-alternatives --install "/usr/bin/gcc" "gcc" "$(which gcc-10)" 100 \
                         --slave "/usr/bin/g++" "g++" "$(which g++-10)" \
                         --slave "/usr/bin/gfortran" "gfortran" "$(which gfortran-10)"
```
If you have other versions of GCC installed on your system, you can run `update-alternatives --install` once for each version, and then execute
```sh
sudo update-alternatives --config gcc
```
to select the correct version.

Alternatively, you could use the standard environment variables `CC`, `CXX` and `FC`, but I have not tried this.

## Prepare virtual environment

From the root of the project directory tree:

```sh
python3.9 -m venv py-venv
. ./py-venv/bin/activate
```

## Install Python dependencies

```sh
pip install -r scripts/requirements.txt
```
## Install Eigen, CasADi, Ipopt, CUTEst, yaml-cpp

```sh
./scripts/install-openblas.sh   # https://www.openblas.net/
./scripts/install-eigen.sh      # https://eigen.tuxfamily.org/index.php
./scripts/install-casadi.sh     # https://web.casadi.org/
./scripts/install-cutest.sh     # https://github.com/ralna/CUTEst
./scripts/install-lbfgspp.sh    # https://github.com/yixuan/LBFGSpp
./scripts/install-yaml-cpp.sh   # https://github.com/jbeder/yaml-cpp
./scripts/install-gtest.sh      # https://google.github.io/googletest/
```

## Build tests

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Asan
cmake --build . -j$(nproc) -t tests
./test/tests
```