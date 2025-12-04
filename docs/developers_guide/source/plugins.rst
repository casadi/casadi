Developing Solver Plugins
=========================

.. include:: defs.rst

What follows are some initial notes on how to develop a |casadi| plugin for your solver.

General Notes
-------------

* In |casadi| your solver plugin is essentially a stateful |Function| of your initial point, parameter, and bounds.
  In the case of a `conic` plugin (which encompasses (possibly mixed-integer linear and quadratic solvers).
*

`NlpsolPlugin`
--------------
