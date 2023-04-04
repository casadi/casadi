# Installation instructions {#installation}

This page contains detailed instructions for building and installing all of
alpaqa from scratch.
For the installation of the Python package without building it from source,
please see [these instructions](../Sphinx/install/installation.html).

## Linux

### Tools
First, install some basic tools: C and C++ compilers, Git, and Python
(you'll need the development version to build alpaqa's Python interface, and we
install the `venv` module to create virtual environments).
```sh
sudo apt install g++ gcc git python3-venv python3-dev
```
The alpaqa package requires a relatively recent compiler
(tested using GCC 10-12, Clang 14-16).

To install GCC 11 on older versions of Ubuntu, you can use
```sh
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt install gcc-11 g++-11
```
To install the latest version of Clang, you can use the instructions from <https://apt.llvm.org>:
```sh
bash -c "$(wget -O- https://apt.llvm.org/llvm.sh)"
```

### Clone the repository

```sh
git clone https://github.com/kul-optec/alpaqa
```

### Create a virtual environment

For convenience, we'll install everything into a Python virtual environment
(including the C++ libraries and dependencies). This allows you to easily
experiment in a sandbox, without requiring root permissions, and without the
risk of messing with system packages.

```sh
cd alpaqa
python3 -m venv py-venv
. ./py-venv/bin/activate
pip install cmake ninja casadi numpy
```

### Install dependencies

The `scripts` folder contains some Bash scripts to install the necessary
dependencies. Feel free to inspect and modify the installation scripts.
If you already have the dependencies installed globally you can skip these
steps.

```sh
bash ./scripts/install-casadi-static.sh "$VIRTUAL_ENV" Release
bash ./scripts/install-gtest.sh "$VIRTUAL_ENV" Release
bash ./scripts/install-eigen.sh "$VIRTUAL_ENV" Release
```

CasADi is built as a static library because it is later statically linked into
the final alpaqa libraries for better portability, this especially useful when
creating the Python package. If you need to link against CasADi dynamically, you
can use the `install-casadi.sh` script instead.

### Build and install

The following commands build and install the alpaqa C++ library into the virtual
environment.  
You may want to change the installation prefix, e.g. use `--prefix /usr/local`
for a system-wide install (requires `sudo`), or `--prefix $HOME/.local` to
install it for the current user.

```sh
cmake -S. -Bbuild -G "Ninja Multi-Config"

cmake --build build --config Release -j       # Build in release mode
cmake --build build -t test                   # Run the tests
cmake --install build --prefix "$VIRTUAL_ENV" # Install the release version

cmake --build build --config Debug -j         # Build in debug mode
cmake --build build -t test                   # Run the tests with extra checks
cmake --install build --prefix "$VIRTUAL_ENV" # Install the debug version
```
Installing both the release and debug versions can be very useful for checking
matrix dimension errors and out of bounds accesses during development, and 
switching to an optimized version later.

> **Note**  
> If you changed the installation prefix, and
> unless you installed the package to a system folder like `/usr/local`, you'll
> have to add `~/.local` to the `CMAKE_PREFIX_PATH`, e.g. by adding the
> following to your `~/.profile` file, where `$HOME/.local` was the prefix used
> in the when installing alpaqa earlier:
> ```sh
> export CMAKE_PREFIX_PATH="$HOME/.local:$CMAKE_PREFIX_PATH"
> ```
> Then source it (`. ~/.profile`) or log out and back in again.

## Windows

==TODO==

The instructions for Windows are quite similar to the ones for Linux. To install
the dependencies, you can use the Powershell scripts instead of the Bash scripts:

```ps1
./scripts/install-casadi-static.ps1
./scripts/install-gtest.ps1
./scripts/install-eigen.ps1
```

## macOS

==TODO==

The instructions for macOS are the same as the ones for Linux, with the caveat
that the default AppleClang compiler might not be supported. Instead, it is
recommended to use a mainline Clang compiler (version 14 or higher).  
You can select the compiler to use by setting the `CC` and `CXX` environment
variables, for example:
```sh
export CC=clang-14
export CXX=clang++-14
```

***

# Usage

Once the library is installed, you can use it in your own projects.

For example:

**main.cpp**
```cpp
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>

int main() {
    // Use the solvers as shown in the examples
}
```

**CMakeLists.txt**
```cmake
cmake_minimum_required(VERSION 3.16)
project(Project)

# Find the library you just installed:
find_package(alpaqa REQUIRED)

add_executable(main main.cpp)
# Link your executable with the library:
target_link_libraries(main PRIVATE alpaqa::alpaqa)
```

# Python

After creating the virtual environment and installing the dependencies, you can
install the Python module using:
```sh
pip install .
```
To build the Python package without installing, you can use:
```sh
pip install build
python3 -m build .
```
