# Examples {#examples}

The following is a list of examples demonstrating how to use the alpaqa library. More Python examples can be found in the [Sphinx documentation](../Sphinx/examples/index.html). An example and tips on how to structure a C++ project that uses alpaqa can be found in [examples/CMake/Solver](https://github.com/kul-optec/alpaqa/tree/develop/examples/CMake/Solver).

## C++ examples {#cpp_examples}

- @ref C++/CustomCppProblem/main.cpp : How to define and solve a custom problem in C++
- @ref C++/SimpleUnconstrProblem/main.cpp : Similar, but for problems without constraints
- @ref C++/FortranProblem/main.cpp : How to define a custom problem in Fortran, and how to solve it in C++
- @ref C++/DLProblem/main.cpp : How to define a problem that can be loaded dynamically, and how to solve it in C++
- @ref C++/CasADi/Rosenbrock/main.cpp : How to load a pre-compiled CasADi problem, and how to solve it in C++
- @ref C++/CustomControlCppProblem/main.cpp : How to define and solve a custom optimal control problem (OCP) in C++
- @ref problems/sparse-logistic-regression.cpp : How to define a problem that can be loaded dynamically, and how to solve it using the `alpaqa-driver` program
- @ref C++/Advanced/lasso-fbs.cpp : How to use the @ref alpaqa::prox functions

## CMake examples {#cmake_examples}

- [**CMake/Solver**](https://github.com/kul-optec/alpaqa/tree/develop/examples/CMake/Solver) : How te set up a Conan + CMake project that uses alpaqa

## Python examples {#python_examples}

- @ref Python/simple_optimization/getting-started.py : Define and solve a simple nonlinear program in Python
- @ref Python/simple_optimization/rosenbrock.py : Similar, but also includes a visualization of the iterates
- More examples in the [**Sphinx documentation**](../Sphinx/examples/index.html) ...

## Matlab examples {#matlab_examples}

- [**Matlab/getting_started.m**](../Sphinx/reference/matlab-api.html)<!-- not rendered correctly without this comment on Doxygen -->
