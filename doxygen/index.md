# alpaqa

`alpaqa` is an efficient implementation of an Augmented Lagrangian method for general nonlinear programming problems.
It makes use of the first-order, matrix-free PANOC algorithm as an inner solver.
The numerical algorithms themselves are implemented in C++ for optimal 
performance, and they are also exposed as an easy-to-use Python package.

The solvers in this library solve minimization problems of the following form:
@f[
\begin{aligned}
    & \underset{x}{\text{minimize}}
    & & f(x) &&&& f : \Rn \rightarrow \R \\
    & \text{subject to}
    & & \underline{x} \le x \le \overline{x} \\
    &&& \underline{z} \le g(x) \le \overline{z} &&&& g : \Rn \rightarrow \R^m
\end{aligned}
@f]

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