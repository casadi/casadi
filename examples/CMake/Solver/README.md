# alpaqa CMake example: Solver usage

This is a CMake example project that makes use of the alpaqa solvers. It
demonstrates how to easily build and install all dependencies using Conan, and
how you can add the alpaqa libraries to your CMake configuration.

The project itself is very simple:
1. A problem class is defined in `problem.{hpp,cpp}`, which represents a
   standard quadratic program (QP) formulation.
2. The matrices and vectors comprising the QP are loaded from the CSV files in
   the `data` directory.
3. The main function defined in `alpaqa-qp-solver.cpp` includes one of the
   alpaqa solvers, configures some parameters, and then applies it to the QP.

## Before you start

You'll need a recent C++ compiler (e.g. GCC >= 11, Clang >= 16, MSVC >= 19.30)
and Conan. The latter can be installed using Pip:

```sh
pip install -U conan
```

## Instructions

Export the alpaqa library to Conan:

```sh
git clone https://github.com/kul-optec/alpaqa --branch=1.0.0a19 --single-branch
conan export alpaqa
```

Download Conan recipes for the dependencies:

```sh
git clone https://github.com/tttapa/conan-recipes
conan remote add tttapa-conan-recipes "$PWD/conan-recipes"
```

Build and install the dependencies for the example project:

```sh
cd alpaqa/examples/CMake/Solver
conan install . --build=missing -c tools.build:skip_test=True
```

Configure and build the example project:
```sh
cmake --preset conan-release  # Linux, macOS (single-configuration generator)
cmake --preset conan-default  # Windows (multi-configuration generator)
cmake --build --preset conan-release
```

Run the example:

```sh
./build/Release/alpaqa-qp-solver data
```

---

## Advanced

### Shared libraries

To link to the alpaqa shared libraries instead of statically, use the following
Conan options:

```sh
conan install . --build=missing -o 'alpaqa/*:shared=True'
```

Remove the CMake cache to ensure that this new version of alpaqa is picked up.

```sh
rm build/Release/CMakeCache.txt  # single-configuration generator
rm build/CMakeCache.txt  # multi-configuration generator
```

### ABI issues

The advantage of Conan is that it will consistently use the same options and the
same versions of dependencies (e.g. Eigen) when building alpaqa itself and when
building your project.

If you decide not to use Conan, it is crucial that you build both alpaqa and
your project with the same version of Eigen, and with the same compiler flags.
If you don't, you will likely end up with ABI incompatibilities, and your code
will not work.

For example, changing the architecture flags (e.g. `-march=skylake`) can change
the Eigen ABI, because different CPUs with different vector extensions require
different alignments.

### CasADi (and other optional components)

If you need optional alpaqa features, e.g. the CasADi interface, you'll have to
specify them explicitly.

For example, to enable CasADi, first export CasADi itself, and then rebuild
alpaqa with CasADi support enabled:
```sh
conan install . --build=missing -o 'alpaqa/*:with_external_casadi=True'
```

If your project always requires CasADi, add the following options at the bottom
of your project's `conanfile.txt`:
```conanfile
[options]
alpaqa/*:with_external_casadi=True
```

In your project's `CMakeLists.txt`, enable the optional CasADi component:

```cmake
# Locate the alpaqa package
find_package(alpaqa 1.0.0 REQUIRED COMPONENTS CasADi)
# Alternatively, use OPTIONAL_COMPONENTS if CasADi is optional
```

You can now use the alpaqa CasADi loader target:
```cmake
target_link_libraries(alpaqa-qp-solver PRIVATE alpaqa::casadi-loader)
```

A complete list of the available components and targets can be found on the
[CMake API Reference](https://kul-optec.github.io/alpaqa/1.0.0a19/Sphinx/reference/cmake-api.html) page.
