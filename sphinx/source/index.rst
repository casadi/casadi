Alpaqa
===================================
.. include:: shared/definitions.rst

|pylib_name| is an efficient implementation of an Augmented Lagrangian solver for general nonlinear programming problems,
which uses the first-order, matrix-free PANOC algorithm as an inner solver.  
The numerical algorithms themselves are implemented in C++ for optimal performance,
and they are exposed as an easy-to-use Python package.

The solvers in this library solve minimization problems of the following form:

.. math::
    \newcommand\mymathbb[1]
    { {\rm I\mathchoice{\hspace{-2pt}}{\hspace{-2pt}}
        {\hspace{-1.75pt}}{\hspace{-1.7pt}}#1} }
    \newcommand{\Re}{\mymathbb R}
    \begin{aligned}
        & \underset{x}{\text{minimize}}
        & & f(x) &&&& f : \Re^n \rightarrow \Re \\
        & \text{subject to}
        & & \underline{x} \le \phantom{g(}x\phantom{)} \le \overline{x} \\
        &&& \underline{z} \le g(x) \le \overline{z} &&&& g : \Re^n \rightarrow \Re^m
    \end{aligned}

The objective function :math:`f(x)` and the constraints function :math:`g(x)`
should have a Lipschitz-continuous gradient.

For more information, please see :ref:`getting started`.

.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:

   install/*
   usage/*
   examples/examples_landing_page
   reference/python-api.rst
   reference/cpp-api.rst
   contribution/*
   changelog