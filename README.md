[![CI Tests: C++](https://github.com/kul-optec/alpaqa/actions/workflows/linux.yml/badge.svg)](https://github.com/kul-optec/alpaqa/actions/workflows/linux.yml)
[![CI Tests: Python](https://github.com/kul-optec/alpaqa/actions/workflows/wheel-short-test.yml/badge.svg)](https://github.com/kul-optec/alpaqa/actions/workflows/wheel-short-test.yml)
[![CI: Matlab](https://github.com/kul-optec/alpaqa/actions/workflows/matlab.yml/badge.svg)](https://github.com/kul-optec/alpaqa/actions/workflows/matlab.yml)
[![Docs](https://img.shields.io/badge/Documentation-1.1.0a1-blue?logo=sphinx)](https://kul-optec.github.io/alpaqa/1.1.0a1/Sphinx/index.html)
[![PyPI Downloads](https://img.shields.io/pypi/dm/alpaqa?label=PyPI&logo=python)](https://pypi.org/project/alpaqa)
[![GitHub all releases](https://img.shields.io/github/downloads/kul-optec/alpaqa/total?label=Downloads&logo=cplusplus)](https://github.com/kul-optec/alpaqa/releases)
[![GitHub License](https://img.shields.io/github/license/kul-optec/alpaqa?label=License&logo=gnu)](https://github.com/kul-optec/alpaqa/blob/develop/LICENSE)
<!-- ![GitHub top language](https://img.shields.io/github/languages/top/kul-optec/alpaqa) -->

# alpaqa

`alpaqa` is an efficient implementation of an augmented Lagrangian method for
general nonlinear programming problems, which uses the first-order, matrix-free
PANOC algorithm as an inner solver.  
The numerical algorithms themselves are implemented in C++ for optimal
performance, and they are exposed as an easy-to-use Python package. An
experimental MATLAB interface is available as well.

The solvers in this library solve minimization problems of the following form:

$$
    \begin{equation}
        \begin{aligned}
            & \underset{x}{\textbf{minimize}}
            & & f(x) &&&& f : {\rm I\\!R}^n \rightarrow {\rm I\\!R} \\
            & \textbf{subject to}
            & & \underline{x} \le x \le \overline{x} \\
            &&& \underline{z} \le g(x) \le \overline{z} &&&& g : {\rm I\\!R}^n \rightarrow {\rm I\\!R}^m
        \end{aligned}
    \end{equation}
$$

## Documentation

- [**Documentation** (Sphinx)](https://kul-optec.github.io/alpaqa/1.1.0a1/Sphinx/index.html)
- [**Python examples**](https://kul-optec.github.io/alpaqa/1.1.0a1/Sphinx/examples/index.html)
- [**C++ documentation** (Doxygen)](https://kul-optec.github.io/alpaqa/1.1.0a1/Doxygen/index.html)
- [**C++ examples**](https://kul-optec.github.io/alpaqa/1.1.0a1/Doxygen/examples.html)
- [**Matlab documentation**](https://kul-optec.github.io/alpaqa/1.1.0a1/Sphinx/reference/matlab-api.html)

## Installation

The Python interface can be installed directly from [PyPI](https://pypi.org/project/alpaqa):

```sh
python3 -m pip install --upgrade --pre alpaqa
```

For more information, please see the full
[installation instructions](https://kul-optec.github.io/alpaqa/1.1.0a1/Sphinx/install/installation.html).

## Publications

> Pieter Pas, Mathijs Schuurmans, and Panagiotis Patrinos. [Alpaqa: A matrix-free solver for nonlinear MPC and large-scale nonconvex optimization](https://ieeexplore.ieee.org/document/9838172/). In _2022 European Control Conference (ECC)_, pages 417–422, 2022. 
