#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2024 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
"""Coverage for the numpy ufunc protocol (NEP 13) and array-function
protocol (NEP 18) bridges in casadi.  Each test was chosen to lock in a
case the previous string-mangling dispatcher silently mishandled."""

import unittest
import numpy as np
import casadi
from casadi import SX, MX, DM, Function, vertcat
from helpers import casadiTestCase


def _evalf_to_dm(expr):
    """Evaluate a closed-form expression (no symbols) to a DM."""
    if isinstance(expr, DM):
        return expr
    if isinstance(expr, (int, float)):
        return DM(expr)
    return Function("f", [], [expr])()['o0']


class NumpyUfuncTests(casadiTestCase):
    """NEP-13 __array_ufunc__ dispatch.

    Each unary/binary check picks the canonical numpy ufunc and asserts
    a casadi value comes back, matching the analytic answer.  Many of
    these previously raised RuntimeWarning("Implicit conversion of
    symbolic CasADi type to numeric matrix not supported") on the old
    string-mangling dispatcher."""

    # ------------------------- unary ufuncs ------------------------- #

    def _check_unary(self, ufunc, casadi_fn, x_val):
        """Check ufunc(SX)/ufunc(MX)/ufunc(DM) against casadi_fn(x_val)."""
        for cls in (SX, MX):
            x = cls.sym("x")
            f = Function("f", [x], [ufunc(x)])
            ref = Function("f", [x], [casadi_fn(x)])
            self.checkfunction(f, ref, inputs=[x_val], digits=12)
        self.checkarray(ufunc(DM(x_val)), casadi_fn(DM(x_val)), digits=12)

    def test_unary_negative(self):
        self._check_unary(np.negative, lambda x: -x, 1.7)

    def test_unary_positive(self):
        self._check_unary(np.positive, lambda x: x, 1.7)

    def test_unary_absolute(self):
        self._check_unary(np.absolute, casadi.fabs, -1.7)
        self._check_unary(np.fabs, casadi.fabs, -1.7)

    def test_unary_square(self):
        self._check_unary(np.square, lambda x: x * x, 1.7)

    def test_unary_reciprocal(self):
        self._check_unary(np.reciprocal, lambda x: 1.0 / x, 1.7)

    def test_unary_sqrt(self):
        self._check_unary(np.sqrt, casadi.sqrt, 1.7)

    def test_unary_sign(self):
        self._check_unary(np.sign, casadi.sign, -1.7)

    def test_unary_exp_log(self):
        self._check_unary(np.exp, casadi.exp, 0.5)
        self._check_unary(np.log, casadi.log, 1.7)
        self._check_unary(np.log10, casadi.log10, 1.7)
        self._check_unary(np.log1p, casadi.log1p, 0.7)
        self._check_unary(np.expm1, casadi.expm1, 0.5)

    def test_unary_trig(self):
        for u, c in [(np.sin, casadi.sin), (np.cos, casadi.cos),
                     (np.tan, casadi.tan), (np.arctan, casadi.atan)]:
            self._check_unary(u, c, 0.7)
        self._check_unary(np.arcsin, casadi.asin, 0.5)
        self._check_unary(np.arccos, casadi.acos, 0.5)

    def test_unary_hyperbolic(self):
        for u, c in [(np.sinh, casadi.sinh), (np.cosh, casadi.cosh),
                     (np.tanh, casadi.tanh), (np.arcsinh, casadi.asinh),
                     (np.arctanh, casadi.atanh)]:
            self._check_unary(u, c, 0.5)
        self._check_unary(np.arccosh, casadi.acosh, 1.5)

    def test_unary_floor_ceil(self):
        self._check_unary(np.floor, casadi.floor, 1.7)
        self._check_unary(np.ceil, casadi.ceil, 1.3)

    def test_unary_logical_not(self):
        # logic_not is only defined on DM/SX, so just test those.
        self.checkarray(np.logical_not(DM([0, 1, 0])), DM([1, 0, 1]))
        x = SX.sym("x")
        f = Function("f", [x], [np.logical_not(x)])
        self.checkarray(f(0), DM(1))
        self.checkarray(f(1), DM(0))

    # ------------------------- binary ufuncs ------------------------ #

    def _check_binary(self, ufunc, casadi_fn, a_val, b_val):
        for cls in (SX, MX):
            a, b = cls.sym("a"), cls.sym("b")
            f = Function("f", [a, b], [ufunc(a, b)])
            ref = Function("f", [a, b], [casadi_fn(a, b)])
            self.checkfunction(f, ref, inputs=[a_val, b_val], digits=12)
        self.checkarray(ufunc(DM(a_val), DM(b_val)),
                        casadi_fn(DM(a_val), DM(b_val)), digits=12)

    def test_binary_arithmetic(self):
        self._check_binary(np.add,         lambda a, b: a + b, 1.7, 2.3)
        self._check_binary(np.subtract,    lambda a, b: a - b, 1.7, 2.3)
        self._check_binary(np.multiply,    lambda a, b: a * b, 1.7, 2.3)
        self._check_binary(np.divide,      lambda a, b: a / b, 1.7, 2.3)
        self._check_binary(np.true_divide, lambda a, b: a / b, 1.7, 2.3)
        self._check_binary(np.power,       lambda a, b: a ** b, 1.7, 2.3)

    def test_binary_minmax(self):
        for u, c in [(np.fmin,    casadi.fmin), (np.fmax, casadi.fmax),
                     (np.minimum, casadi.fmin), (np.maximum, casadi.fmax)]:
            self._check_binary(u, c, 1.7, 2.3)

    def test_binary_misc(self):
        self._check_binary(np.arctan2,   casadi.atan2,    1.7, 2.3)
        self._check_binary(np.hypot,     casadi.hypot,    3.0, 4.0)
        self._check_binary(np.copysign,  casadi.copysign, 1.7, -1.0)

    def test_binary_compare(self):
        for u, c in [(np.less,          lambda a, b: a < b),
                     (np.less_equal,    lambda a, b: a <= b),
                     (np.greater,       lambda a, b: a > b),
                     (np.greater_equal, lambda a, b: a >= b),
                     (np.equal,         lambda a, b: a == b),
                     (np.not_equal,     lambda a, b: a != b)]:
            self._check_binary(u, c, 1.0, 2.0)

    def test_mod_numpy_semantics(self):
        # numpy.mod / numpy.remainder follow Python %: sign of divisor.
        # casadi.fmod and casadi.remainder do NOT.  The dispatch must
        # bridge to the right semantics.  Both np.mod and np.remainder
        # are the same ufunc whose __name__ is "remainder".
        cases = [(-3, 5, 2), (-7, 5, 3), (-3, -5, -3), (3, -5, -2),
                 (7, -3, -2), (-3.5, 2, 0.5)]
        for a, b, expected in cases:
            self.checkarray(np.mod(DM(a), DM(b)),       DM(expected), digits=12)
            self.checkarray(np.remainder(DM(a), DM(b)), DM(expected), digits=12)
        # symbolic
        a, b = SX.sym("a"), SX.sym("b")
        f = Function("f", [a, b], [np.mod(a, b)])
        for a_val, b_val, expected in cases:
            self.checkarray(f(a_val, b_val), DM(expected), digits=12)

    def test_fmod_c_semantics(self):
        # numpy.fmod follows C: sign of dividend.
        cases = [(-3, 5, -3), (-7, 5, -2), (3, -5, 3), (-3, -5, -3),
                 (7, -3, 1), (-3.5, 2, -1.5)]
        for a, b, expected in cases:
            self.checkarray(np.fmod(DM(a), DM(b)), DM(expected), digits=12)

    def test_logical_ops(self):
        a = DM([1, 0, 1, 0])
        b = DM([1, 1, 0, 0])
        self.checkarray(np.logical_and(a, b), DM([1, 0, 0, 0]))
        self.checkarray(np.logical_or(a, b),  DM([1, 1, 1, 0]))
        self.checkarray(np.logical_not(a),    DM([0, 1, 0, 1]))

    # --------------------- reduce / accumulate ---------------------- #

    def test_sum_reduce(self):
        for cls in (SX, MX):
            x = cls.sym("x", 4)
            f = Function("f", [x], [np.sum(x)])
            ref = Function("f", [x], [casadi.sum(x)])
            self.checkfunction(f, ref, inputs=[[1.0, 2.0, 3.0, 4.0]], digits=14)
        # DM stays DM (not float) on the new dispatcher
        self.assertIsInstance(np.sum(DM([1, 2, 3])), DM)
        self.checkarray(np.sum(DM([1, 2, 3])), DM(6))
        self.checkarray(np.sum(DM([[1, 2, 3], [4, 5, 6]])), DM(21))

    def test_sum_axis(self):
        x = DM([[1, 2, 3], [4, 5, 6]])
        # axis=0: sum down columns -> 1xN
        self.checkarray(np.sum(x, axis=0), DM([[5, 7, 9]]))
        # axis=1: sum across rows -> Nx1
        self.checkarray(np.sum(x, axis=1), DM([6, 15]))

    def test_cumsum_accumulate(self):
        for cls in (SX, MX):
            x = cls.sym("x", 4)
            f = Function("f", [x], [np.cumsum(x)])
            g = Function("g", [x], [np.add.accumulate(x)])
            self.checkfunction(f, g, inputs=[[1.0, 2.0, 3.0, 4.0]], digits=14)
        # numeric
        self.checkarray(np.cumsum(DM([1, 2, 3])), DM([1, 3, 6]))

    def test_minmax_reduce(self):
        x = DM([[1, 5, 3], [2, 4, 6]])
        self.checkarray(np.max(x), DM(6))
        self.checkarray(np.min(x), DM(1))
        self.checkarray(np.maximum.reduce(x.nonzeros()), DM(6))
        # global only -- per-axis max/min should fall through cleanly
        with self.assertRaises(Exception):
            np.max(x, axis=0)

    def test_all_any_reduce(self):
        self.checkarray(np.all(DM([1, 1, 1])), DM(1))
        self.checkarray(np.all(DM([1, 0, 1])), DM(0))
        self.checkarray(np.any(DM([0, 0, 0])), DM(0))
        self.checkarray(np.any(DM([0, 0, 1])), DM(1))

    # ----------------------------- mixed --------------------------- #

    def test_mixed_numpy_array(self):
        # Plain numpy arrays mixing with casadi types route through the
        # same ufunc dispatch.
        x = SX.sym("x")
        result = np.add(np.array([1.0, 2.0]), x)
        f = Function("f", [x], [result])
        self.checkarray(f(3.0), DM([4, 5]).T if False else DM([4, 5]))

    def test_no_silent_warnings(self):
        # The old dispatcher emitted RuntimeWarning("Implicit conversion
        # of symbolic CasADi type to numeric matrix not supported") for
        # every unrecognized ufunc.  Make sure we don't.
        import warnings
        x = SX.sym("x")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            for op in (np.negative, np.square, np.reciprocal,
                       np.absolute, np.expm1, np.log1p):
                op(x)


class NumpyArrayFunctionTests(casadiTestCase):
    """NEP-18 __array_function__ dispatch.  Lets users write idiomatic
    numpy code that just-works on casadi types."""

    # ------------------------- stacking ----------------------------- #

    def test_concatenate(self):
        a = DM([[1, 2], [3, 4]])
        b = DM([[5, 6], [7, 8]])
        self.checkarray(np.concatenate([a, b], axis=0),
                        DM([[1, 2], [3, 4], [5, 6], [7, 8]]))
        self.checkarray(np.concatenate([a, b], axis=1),
                        DM([[1, 2, 5, 6], [3, 4, 7, 8]]))

    def test_vstack_hstack(self):
        a, b = DM([[1, 2]]), DM([[3, 4]])
        self.checkarray(np.vstack([a, b]), DM([[1, 2], [3, 4]]))
        self.checkarray(np.hstack([a, b]), DM([[1, 2, 3, 4]]))

    # ------------------------- shape -------------------------------- #

    def test_reshape_C_order(self):
        # numpy's default order='C' (row-major).  casadi.reshape is
        # F-order, so the bridge transposes to get numpy semantics.
        v = DM([1, 2, 3, 4, 5, 6])
        self.checkarray(np.reshape(v, (2, 3)), DM([[1, 2, 3], [4, 5, 6]]))
        self.checkarray(np.reshape(v, (3, 2)), DM([[1, 2], [3, 4], [5, 6]]))
        # Explicit F-order matches casadi.reshape
        self.checkarray(np.reshape(v, (2, 3), order='F'),
                        DM([[1, 3, 5], [2, 4, 6]]))

    def test_transpose(self):
        m = DM([[1, 2, 3], [4, 5, 6]])
        self.checkarray(np.transpose(m), DM([[1, 4], [2, 5], [3, 6]]))

    # ------------------------- linear algebra ----------------------- #

    def test_dot_matmul(self):
        # vector dot
        self.checkarray(np.dot(DM([1, 2, 3]), DM([4, 5, 6])), DM(32))
        # matrix-matrix
        A = DM([[1, 2], [3, 4]])
        B = DM([[5, 6], [7, 8]])
        ref = DM([[19, 22], [43, 50]])
        self.checkarray(np.dot(A, B), ref)
        self.checkarray(np.matmul(A, B), ref)
        # symbolic matmul (np.matmul is a ufunc, dispatched via NEP 13)
        for cls in (SX, MX):
            X, Y = cls.sym("X", 2, 2), cls.sym("Y", 2, 2)
            f = Function("f", [X, Y], [np.matmul(X, Y)])
            ref = Function("ref", [X, Y], [casadi.mtimes(X, Y)])
            self.checkfunction(f, ref,
                               inputs=[A.full(), B.full()], digits=12)

    def test_inner_outer(self):
        a, b = DM([1, 2, 3]), DM([4, 5, 6])
        self.checkarray(np.inner(a, b), DM(32))
        self.checkarray(np.outer(a, b),
                        DM([[4, 5, 6], [8, 10, 12], [12, 15, 18]]))

    def test_kron_cross(self):
        self.checkarray(np.kron(DM([[1, 2]]), DM([[3], [4]])),
                        DM([[3, 6], [4, 8]]))
        self.checkarray(np.cross(DM([1, 0, 0]), DM([0, 1, 0])),
                        DM([0, 0, 1]))

    def test_diag(self):
        self.checkarray(np.diag(DM([1, 2, 3])), DM([[1, 0, 0], [0, 2, 0], [0, 0, 3]]))

    def test_linalg_norm(self):
        self.checkarray(np.linalg.norm(DM([3, 4])), DM(5))
        self.checkarray(np.linalg.norm(DM([3, -4]), ord=1), DM(7))
        self.checkarray(np.linalg.norm(DM([3, -4]), ord=np.inf), DM(4))

    def test_linalg_solve(self):
        A = DM([[2, 0], [0, 3]])
        b = DM([4, 9])
        self.checkarray(np.linalg.solve(A, b), DM([2, 3]))

    def test_linalg_det_inv(self):
        A = DM([[1, 2], [3, 4]])
        self.checkarray(np.linalg.det(A), DM(-2))
        self.checkarray(np.linalg.inv(A), DM([[-2, 1], [1.5, -0.5]]))

    def test_linalg_cholesky(self):
        A = DM([[4, 2], [2, 5]])
        L = np.linalg.cholesky(A)            # lower-triangular by numpy convention
        self.checkarray(L @ L.T, A, digits=12)

    # ------------------------- aerosandbox-inspired ----------------- #

    def test_clip(self):
        self.checkarray(np.clip(DM([-1, 0.5, 2]), 0, 1), DM([0, 0.5, 1]))
        x = SX.sym("x")
        f = Function("f", [x], [np.clip(x, -1, 1)])
        self.checkarray(f(-2), DM(-1))
        self.checkarray(f(0.3), DM(0.3))
        self.checkarray(f(5), DM(1))

    def test_diff(self):
        self.checkarray(np.diff(DM([1, 3, 6, 10])), DM([2, 3, 4]))
        x = SX.sym("x", 4)
        f = Function("f", [x], [np.diff(x)])
        self.checkarray(f([1, 3, 6, 10]), DM([2, 3, 4]))

    def test_roll(self):
        # 1-D wrap
        self.checkarray(np.roll(DM([1, 2, 3, 4, 5]),  2), DM([4, 5, 1, 2, 3]))
        self.checkarray(np.roll(DM([1, 2, 3, 4, 5]), -1), DM([2, 3, 4, 5, 1]))
        # 2-D axis
        m = DM([[1, 2, 3], [4, 5, 6]])
        self.checkarray(np.roll(m, 1, axis=0), DM([[4, 5, 6], [1, 2, 3]]))
        self.checkarray(np.roll(m, 1, axis=1), DM([[3, 1, 2], [6, 4, 5]]))

    # ------------------------- conditional -------------------------- #

    def test_where(self):
        x = SX.sym("x")
        result = np.where(x > 0, x, -x)
        f = Function("f", [x], [result])
        self.checkarray(f(2.0),  DM(2.0))
        self.checkarray(f(-3.0), DM(3.0))
        # numeric
        self.checkarray(np.where(DM([1, -1, 1]) > 0, DM([10, 20, 30]), DM([-10, -20, -30])),
                        DM([10, -20, 30]))

    # ------------------------- like-constructors -------------------- #

    def test_zeros_ones_like(self):
        a = DM([[1, 2], [3, 4]])
        self.checkarray(np.zeros_like(a), DM([[0, 0], [0, 0]]))
        self.checkarray(np.ones_like(a),  DM([[1, 1], [1, 1]]))
        self.assertIsInstance(np.zeros_like(SX.sym("x", 2, 2)), SX)
        self.assertIsInstance(np.ones_like(MX.sym("x", 2, 2)), MX)

    def test_linspace(self):
        self.checkarray(np.linspace(DM(0), DM(1), 5),
                        DM([0, 0.25, 0.5, 0.75, 1.0]))

    def test_tile(self):
        self.checkarray(np.tile(DM([[1, 2]]), (2, 3)),
                        DM([[1, 2, 1, 2, 1, 2], [1, 2, 1, 2, 1, 2]]))

    # ------------------------- failure modes ------------------------ #

    def test_unsupported_returns_clear_error(self):
        # Symbolic operations that have no casadi equivalent must
        # produce a clear failure (not a "implicit conversion" warning).
        with self.assertRaises(Exception):
            np.all(MX([1, 0, 1]))
        # Unknown numpy functions return NotImplemented from our
        # dispatch and let numpy raise its own clean message.
        with self.assertRaises(Exception):
            np.argmax(SX.sym("x", 3))


if __name__ == '__main__':
    unittest.main()
