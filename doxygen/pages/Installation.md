# Installation instructions {#installation}

> **Note**  
> This page contains detailed instructions for building and installing all of
> alpaqa from source.
> For the installation of the Python package without building it from source,
> as well as installation instructions for pre-built, released binaries of the
> C++ library and the MATLAB interface, please see [these instructions](../Sphinx/install/installation.html)
> instead.
> For instructions on how to get alpaqa through the [Conan](https://conan.io/)
> package manager, see the @ref cmake_examples.

## Linux

### Tools
First, install some basic tools: C and C++ compilers, Git, and Python
(you'll need the development version to build alpaqa's Python interface, and we
install the `venv` module to create virtual environments).
```sh
sudo apt install build-essential git python3-venv python3-dev
```
The alpaqa package requires a relatively recent compiler
(tested using GCC 11-15, Clang (libc++) 16-21, or Clang (libstdc++) 17-21).

### Clone the repository

```sh
git clone https://github.com/kul-optec/alpaqa --branch=develop --single-branch
```

### Create a virtual environment

For convenience, we'll install everything into a Python virtual environment.
This allows you to easily experiment in a sandbox, without requiring root
permissions, and without the risk of messing with system packages.

```sh
cd alpaqa
python3 -m venv .venv
. ./.venv/bin/activate
pip install conan
```

If you haven't used Conan before, run the following command to configure a
default Conan profile for your system:
```sh
conan profile detect
```

### Install dependencies

The Conan package manager is used for installing the dependencies. Since not all
packages are in the main Conan Center repository, we add a secondary repository:

```sh
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes" --force
```

We'll now install the dependencies of alpaqa. See `conanfile.py` for a list of
the available options. For the sake of simplicity, we'll simply install the
dependencies for the default, minimal configuration of alpaqa:

```sh
conan install . --build=missing -s build_type=Release \
    -c tools.cmake.cmaketoolchain:generator="Ninja Multi-Config" \
    -c tools.build:skip_test=True
```

### Build and install

The following commands build and install the alpaqa C++ library into the virtual
environment.  
You may want to change the installation prefix, e.g. use `--prefix /usr/local`
for a system-wide install (requires `sudo`), or `--prefix $HOME/.local` to
install it for the current user.

```sh
cmake --preset conan-default
cmake --build --preset conan-release
cmake --install build --prefix "$VIRTUAL_ENV" --config Release  # optional
```

Alternatively, avoid manual installation by consuming alpaqa through Conan.
To make the `alpaqa` package available, use:
```sh
conan export .
```
See the @ref cmake_examples "`examples/CMake/Solver` example project"
for details.

## Windows

The instructions for Windows are quite similar to the ones for Linux, but you'll
usually want to use the default generator rather than Ninja when doing
`conan install`:
```ps1
cd alpaqa
py -m venv .venv
&./.venv/Scripts/Activate.ps1
pip install conan
conan profile detect
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes" --force
conan install . --build=missing -s build_type=Release -c tools.build:skip_test=True
cmake --preset conan-default
cmake --build --preset conan-release
cmake --install build --prefix "$env:VIRTUAL_ENV" --config Release  # optional
```

## macOS

The instructions for macOS are the same as the ones for Linux, with the caveat
that the default AppleClang compiler might not yet support the necessary C++20
features used by alpaqa. If this is the case, you can use a mainline Clang
compiler (version 16 or higher), that you install using Homebrew or another
package manager.
Make sure you use a Conan profile that selects the appropriate version of Clang.

***

# Usage

Once the library is installed, you can use it in your own projects.

It is highly recommended to have a look at the
@ref cmake_examples "`examples/CMake/Solver` example project",
but in short:

**main.cpp**
```cpp
#include <alpaqa/panoc-alm.hpp>

int main() {
    // Use the solvers as shown in the examples
}
```

**CMakeLists.txt**
```cmake
cmake_minimum_required(VERSION 3.17...4.1)
project(Project)

# Find the library you just installed:
find_package(alpaqa REQUIRED)

add_executable(main main.cpp)
# Link your executable with the library:
target_link_libraries(main PRIVATE alpaqa::alpaqa)
```

Different components and targets are available. Depending on your needs, you
might want to link to:

 - `alpaqa::alpaqa`: the core alpaqa library and solvers
 - `alpaqa::casadi-loader`: provides the `CasADiProblem` class that allows
    the solvers to interface with problems formulated using CasADi
 - `alpaqa::casadi-ocp-loader`: experimental optimal-control specific CasADi
    problem specification
 - `alpaqa::dl-api`: the stand-alone C API for formulating problems that can be
    loaded dynamically by alpaqa (`alpaqa/dl/dl-problem.h`)
 - `alpaqa::dl-loader`: provides the `DLProblem` class to load such problems
 - `alpaqa::cutest-interface`: provides the `CUTEstProblem` class for loading
    problems formulated using SIF/CUTEst
 - `alpaqa::ipopt-adapter`: allows passing any alpaqa problem to the Ipopt
    solver
 - `alpaqa::lbfgsb-adapter`: allows passing any alpaqa problem to the L-BFGS-B
    solver
 - `alpaqa::qpalm-adapter`: allows passing any alpaqa problem to the QPALM
    solver

See the [CMake API documentation](../Sphinx/reference/cmake-api.html) for more
details.

# Python

Similarly to the instructions above, create a virtual environment and install
the dependencies using Conan. You want to add the `with_python` option.

## Linux and macOS
```sh
cd alpaqa
python3 -m venv .venv
. ./.venv/bin/activate
pip install conan
conan profile detect ||:
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes" --force
```

## Windows
```ps1
cd alpaqa
py -m venv .venv
&./.venv/Scripts/Activate.ps1
pip install conan
conan profile detect || $null
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes" --force
```

After creating the virtual environment and installing the dependencies, you can
install the Python module using Pip (this may take a while):
```sh
pip install -v ".[test,casadi]"
```
To build the Python package without installing, you can use:
```sh
pip install build
python -m build -w .
```
Finally, test the package:
```sh
pytest
```

# Matlab

## Linux and macOS

```sh
cd alpaqa
python3 -m venv .venv
. ./.venv/bin/activate
python -m pip install -U conan
conan profile detect ||:
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes" --force
conan install . --build=missing \
    -c tools.cmake.cmaketoolchain:generator="Ninja Multi-Config" \
    -c tools.build:skip_test=True \
    -o with_matlab=True -o with_external_casadi=True
cmake --preset conan-matlab
cmake --build --preset conan-matlab-release -t alpaqa_mex
cmake --install build/matlab --config Release \
    --prefix ~/Documents/MATLAB --component mex_interface
```

## Windows

```sh
cd alpaqa
py -m venv .venv
&./.venv/Scripts/Activate.ps1
pip install conan
conan profile detect || $null
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes" --force
conan install . --build=missing `
    -c tools.build:skip_test=True `
    -o with_matlab=True -o with_external_casadi=True
cmake --preset conan-matlab
cmake --build --preset conan-matlab-release -t alpaqa_mex
cmake --install build/matlab --config Release `
    --prefix "$env:USERPROFILE\Documents\MATLAB" --component mex_interface
```

## Uninstall

To uninstall the alpaqa MATLAB/MEX interface, simply remove the `+alpaqa`
directory, e.g. by running the following command in the MATLAB command window:

```matlab
rmdir(fullfile(userpath, '+alpaqa'), 's')
```
