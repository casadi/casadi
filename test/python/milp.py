import unittest
from enum import Enum, auto
from typing import Dict, Union

import casadi as ca
from numpy import inf


class AutoName(Enum):
    """Workaround, StrEnum doesn't exist yet till Python 3.11"""

    def _generate_next_value_(name, start, count, last_values):
        return name.lower()


class IntegerSolver(AutoName):
    CBC = auto()
    HIGHS = auto()


def get_solver_options(
    solver: IntegerSolver,
    time_limit: Union[int, float] = None,
) -> Dict[str, Dict]:
    if time_limit is None:
        time_limit = inf
    if solver == IntegerSolver.CBC:
        return {
            solver.value: {
                "loglevel": 0,
                "AllowableGap": 0.0,
                "MaximumSeconds": time_limit,
            }
        }
    elif solver == IntegerSolver.HIGHS:
        return {solver.value: {"output_flag": False, "time_limit": time_limit}}
    else:
        raise ValueError(f"Unknown solver: {solver}")


class TestSolvingMILP(unittest.TestCase):
    """
    Solve simple Mixed-Integer Linear Programming (MILP) problem
    maximize     x + y
    subject to  -2x + 2y >= 1
                -8x + 10y <= 13
                 x,y >= 0
                 x,y are integer

    Optimal solution of LP relaxation is (x, y) = (4, 4.5) with objective 8.5
    Optimal integer solution is (x, y) = (1, 2) with objective 3.0
    """

    def setUp(self):
        # Declare model variables
        x = ca.SX.sym("x")
        y = ca.SX.sym("y")

        # Formulate MILP problem
        self.milp = {
            # Objective (swap sign so we maximize x + y)
            "f": -(x + y),
            # Constraints
            "g": ca.vertcat(-2 * x + 2 * y, -8 * x + 10 * y),
            # Variables
            "x": ca.vertcat(x, y),
        }

        # Define right-hand side values
        self.rhs = {
            "lbg": ca.vertcat(1, -inf),
            "ubg": ca.vertcat(inf, 13),
            "lbx": ca.vertcat(0, 0),
            "ubx": ca.vertcat(inf, inf),
        }

    def test_relaxed_problem(self):
        """Solve the relaxed problem using each solver"""
        for solver in IntegerSolver:
            with self.subTest():
                options = get_solver_options(solver)

                solver_obj = ca.qpsol("nlp", solver.value, self.milp, options)
                results = solver_obj(**self.rhs)

                assert results["f"].full().flatten()[0] == -8.5
                assert (results["x"].full().flatten() == [4, 4.5]).all()

    def test_integer_problem(self):
        """Solve the integer problem using each solver"""
        for solver in IntegerSolver:
            with self.subTest():
                options = get_solver_options(solver)
                options["discrete"] = [True, True]

                solver_obj = ca.qpsol("nlp", solver.value, self.milp, options)
                results = solver_obj(**self.rhs)

                assert results["f"].full().flatten()[0] == -3.0
                assert (results["x"].full().flatten() == [1, 2]).all()

    def test_solver_time_limit(self):
        """Time limit as integer or float should both work"""
        for solver in IntegerSolver:
            for time_limit in [10, 10.0]:
                with self.subTest():
                    options = get_solver_options(solver, time_limit=time_limit)

                    solver_obj = ca.qpsol("nlp", solver.value, self.milp, options)
                    solver_obj(**self.rhs)


if __name__ == "__main__":
    unittest.main()
