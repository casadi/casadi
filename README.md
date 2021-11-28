[![Build Status](https://github.com/kul-optec/alpaqa/workflows/CI%20Tests/badge.svg)](https://github.com/kul-optec/alpaqa/actions)
[![Test Coverage](https://img.shields.io/endpoint?url=https://kul-optec.github.io/alpaqa/Coverage/shield.io.coverage.json)](https://kul-optec.github.io/alpaqa/Coverage/index.html)
[![GitHub](https://img.shields.io/github/stars/kul-optec/alpaqa?label=GitHub&logo=github)](https://github.com/kul-optec/alpaqa)


# alpaqa

`Alpaqa` is an efficient implementation of the Augmented Lagrangian method for general nonlinear programming problems,
which uses the first-order, matrix-free PANOC algorithm as an inner solver.  
The numerical algorithms themselves are implemented in C++ for optimal performance,
and they are exposed as an easy-to-use Python package.

The solvers in this library solve minimization problems of the following form:

<div align="center">

![Problem formulation](https://github.com/kul-optec/alpaqa/blob/main/doxygen/images/problem.svg?raw=True)

</div>

The objective function _f_(x) and the constraints function _g_(x)
should have a Lipschitz-continuous gradient.

## Documentation

[**Sphinx documentation**](https://kul-optec.github.io/alpaqa/Sphinx/index.html)  
[**Doxygen documentation**](https://kul-optec.github.io/alpaqa/Doxygen/index.html)  
[**Python examples**](https://kul-optec.github.io/alpaqa/Sphinx/examples/examples_landing_page.html)  
[**C++ examples**](https://kul-optec.github.io/alpaqa/Doxygen/examples.html)  

## Installation

The project is available on [PyPI](https://pypi.org/project/alpaqa):

```sh
python3 -m pip install alpaqa
```

For more information, please see the full
[installation instructions](https://kul-optec.github.io/alpaqa/Sphinx/install/installation.html).
