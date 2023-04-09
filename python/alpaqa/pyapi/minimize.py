from dataclasses import dataclass
import casadi as cs
from typing import Union, Optional, Tuple
import numpy as np
from copy import copy
from ..casadi_loader import generate_and_compile_casadi_problem
from ..alpaqa import CasADiProblem
import inspect


@dataclass
class MinimizationProblemDescription:
    """
    High-level description of a minimization problem.
    """

    objective_expr: Union[cs.SX, cs.MX]
    variable: Union[cs.SX, cs.MX]
    constraints_expr: Optional[Union[cs.SX, cs.MX]] = None
    penalty_constraints_expr: Optional[Union[cs.SX, cs.MX]] = None
    parameter: Optional[Union[cs.SX, cs.MX]] = None
    parameter_value: Optional[np.ndarray] = None
    regularizer: Optional[Union[float, np.ndarray]] = None
    bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None
    constraints_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None
    penalty_constraints_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None

    @staticmethod
    def _assert_not_set_before(value):
        if value is not None:
            caller_frame = inspect.getouterframes(inspect.currentframe())[1]
            raise ValueError(f"{caller_frame.function} cannot be called twice")

    def subject_to_box(self, C: Tuple[np.ndarray, np.ndarray]):
        """
        Add box constraints :math:`x \\in C` on the problem variables.
        """
        self._assert_not_set_before(self.bounds)
        ret = copy(self)
        ret.bounds = C
        return ret

    def subject_to(
        self,
        g: Union[cs.SX, cs.MX],
        D: Optional[Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]] = None,
    ):
        """
        Add general constraints :math:`g(x) \\in D`, handled using an augmented
        Lagrangian method.
        """
        self._assert_not_set_before(self.constraints_expr)
        self._assert_not_set_before(self.constraints_bounds)
        if D is not None and not isinstance(D, tuple):
            D = D, D
        ret = copy(self)
        ret.constraints_expr = g
        ret.constraints_bounds = D
        return ret

    def subject_to_penalty(
        self,
        g: Union[cs.SX, cs.MX],
        D: Optional[Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]] = None,
    ):
        """
        Add general constraints :math:`g(x) \\in D`, handled using a quadratic
        penalty method.
        """
        self._assert_not_set_before(self.penalty_constraints_expr)
        self._assert_not_set_before(self.penalty_constraints_bounds)
        if D is not None and not isinstance(D, tuple):
            D = D, D
        ret = copy(self)
        ret.penalty_constraints_expr = g
        ret.penalty_constraints_bounds = D
        return ret

    def with_l1_regularizer(self, λ: Union[float, np.ndarray]):
        """
        Add an :math:`\\ell_1`-regularization term :math:`\\|\\lambda x\\|_1`
        to the objective.
        """
        self._assert_not_set_before(self.regularizer)
        ret = copy(self)
        ret.regularizer = λ
        return ret

    def with_param(self, p: Union[cs.SX, cs.MX], value: np.ndarray = None):
        """
        Make the problem depend on a symbolic parameter, with an optional
        default value. The value can be changed after the problem has been
        loaded, as wel as in between solves.
        """
        self._assert_not_set_before(self.parameter)
        ret = copy(self)
        ret.parameter = p
        if value is not None:
            ret.parameter_value = value
        return ret

    def with_param_value(self, value: np.ndarray):
        """
        Explicitly change the parameter value for the parameter added by
        :py:func:`with_param`.
        """
        if self.parameter is None:
            raise RuntimeError("problem has no parameters")
        ret = copy(self)
        ret.parameter_value = value
        return ret

    def compile(self, **kwargs) -> CasADiProblem:
        """
        Generate, compile and load the problem.
        """
        # Function arguments (variables and parameters)
        x, param = self.variable, self.parameter
        args = [cs.vec(x), cs.vec(param)] if param is not None else [cs.vec(x)]
        # Objective and constraints functions
        f = [self.objective_expr]
        g_qpm = self.penalty_constraints_expr
        g_alm = self.constraints_expr
        g = cs.vertcat(*(g for g in (g_qpm, g_alm) if g is not None))
        g = [] if g.shape == (0, 0) else [g]
        # Problem dimensions
        n = cs.vec(x).shape[0]
        p = cs.vec(param).shape[0]
        m_qpm = g_qpm.shape[0] if g_qpm is not None else 0
        m_alm = g_alm.shape[0] if g_alm is not None else 0
        # Bound constraints
        C = self.bounds
        if C is None:
            C = (-np.inf * np.ones(n), +np.inf * np.zeros(n))
        # General quadratic penalty method constraint set
        D_qpm = self.penalty_constraints_bounds
        if D_qpm is None:
            D_qpm = (-np.inf * np.ones(m_qpm), +np.inf * np.zeros(m_qpm))
        # General augmented Lagrangian method constraint set
        D_alm = self.constraints_bounds
        if D_alm is None:
            D_alm = (-np.inf * np.ones(m_alm), +np.inf * np.zeros(m_alm))
        num_param = self.parameter_value
        if num_param is None:
            num_param = np.NaN * np.ones(p)
        λ = self.regularizer
        problem = generate_and_compile_casadi_problem(
            f=cs.Function("f", args, f),
            g=cs.Function("g", args, g),
            C=C,
            D=np.hstack((D_qpm, D_alm)),
            param=num_param,
            l1_reg=λ,
            penalty_alm_split=m_qpm,
            **kwargs,
        )
        return problem


def minimize(
    f: Union[cs.SX, cs.MX], x: Union[cs.SX, cs.MX]
) -> MinimizationProblemDescription:
    """
    Formulate a minimization problem with objective function :math:`f(x)` and
    unknown variables :math:`x`.
    """
    return MinimizationProblemDescription(objective_expr=f, variable=x)


__all__ = [
    "MinimizationProblemDescription",
    "minimize",
]
