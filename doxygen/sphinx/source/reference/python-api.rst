Python API Reference 
====================

High-level problem formulation

.. autofunction:: alpaqa.pyapi.minimize.minimize
    :noindex:

.. autoclass:: alpaqa.pyapi.minimize.MinimizationProblemDescription
    :noindex:
    :no-members:
    :no-special-members:

    .. automethod:: alpaqa.pyapi.minimize.MinimizationProblemDescription.subject_to_box
        :noindex:

    .. automethod:: alpaqa.pyapi.minimize.MinimizationProblemDescription.subject_to
        :noindex:

    .. automethod:: alpaqa.pyapi.minimize.MinimizationProblemDescription.subject_to_penalty
        :noindex:

    .. automethod:: alpaqa.pyapi.minimize.MinimizationProblemDescription.with_l1_regularizer
        :noindex:

    .. automethod:: alpaqa.pyapi.minimize.MinimizationProblemDescription.with_param
        :noindex:

    .. automethod:: alpaqa.pyapi.minimize.MinimizationProblemDescription.with_param_value
        :noindex:

    .. automethod:: alpaqa.pyapi.minimize.MinimizationProblemDescription.compile
        :noindex:

Inner PANOC and ZeroFPR Solvers
-------------------------------

.. autoclass:: alpaqa._alpaqa.float64.PANOCSolver
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.PANOCParams
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.ZeroFPRSolver
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.ZeroFPRParams
    :noindex:

Accelerators
^^^^^^^^^^^^^^^^^^

.. autoclass:: alpaqa._alpaqa.float64.LBFGSDirection
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.StructuredLBFGSDirection
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.LBFGS.Params
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.AndersonDirection
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.AndersonAccel.Params
    :noindex:

Inner PANTR Solver
------------------

.. autoclass:: alpaqa._alpaqa.float64.PANTRSolver
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.PANTRParams
    :noindex:

Accelerators
^^^^^^^^^^^^

.. autoclass:: alpaqa._alpaqa.float64.NewtonTRDirection
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.SteihaugCGParams
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.NewtonTRDirectionParams
    :noindex:

Outer ALM Solver
----------------

.. autoclass:: alpaqa._alpaqa.float64.ALMSolver
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.ALMParams
    :noindex:

Problem formulation
-------------------

.. autoclass:: alpaqa._alpaqa.float64.BoxConstrProblem
    :noindex:
    :inherited-members:

.. autoclass:: alpaqa._alpaqa.float64.Problem
    :noindex:
    :inherited-members:

CasADi Interface
----------------

.. automodule:: alpaqa.casadi_generator
    :noindex:

.. automodule:: alpaqa.casadi_loader
    :noindex:

.. autoclass:: alpaqa._alpaqa.float64.CasADiProblem
    :noindex:

.. autofunction:: alpaqa._alpaqa.float64.load_casadi_problem
    :noindex:

All
---

.. automodule:: alpaqa._alpaqa

.. automodule:: alpaqa._alpaqa.float64

.. automodule:: alpaqa.casadi_generator

.. automodule:: alpaqa.casadi_loader

.. automodule:: alpaqa.pyapi.minimize
