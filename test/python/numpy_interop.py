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
protocol (NEP 18) bridges in casadi.

Every numeric assertion compares the casadi-dispatched call against a
pure numpy reference computed on the same dense input -- never against
a hand-coded literal -- so any drift between casadi semantics and numpy
semantics shows up as a test failure.

Where applicable, each test exercises shape variants (scalar / row
vector / column vector / matrix) and runs against all three casadi
types (DM, SX, MX): the symbolic types build a Function, evaluate it on
the same dense values, and compare to the same numpy reference.  Some
casadi ops have no MX implementation in the core (e.g. det, chol,
logic_all); those tests opt out of MX explicitly via sym_types=(SX,)."""

import unittest
import numpy as np
from casadi import SX, MX, DM, Function, Sparsity
from helpers import casadiTestCase


# ----------------------------- shape fixtures ----------------------------- #

_SCALAR = DM(0.7)
_ROW    = DM([[0.4, 0.6, 0.8]])             # (1, 3)
_COL    = DM([0.4, 0.6, 0.8])               # (3, 1)
_MAT23  = DM([[0.2, 0.4, 0.6], [0.3, 0.5, 0.7]])

_SPD3 = DM([[4.0, 1.0, 0.0],
            [1.0, 3.0, 1.0],
            [0.0, 1.0, 5.0]])
_RHS3 = DM([1.0, 2.0, 3.0])
_RHS3M = DM([[1.0, -2.0],
             [2.0,  0.5],
             [3.0,  1.0]])

_SHAPES_ELEMWISE = [
    ("scalar", _SCALAR),
    ("row3",   _ROW),
    ("col3",   _COL),
    ("mat2x3", _MAT23),
]

# Sparse fixtures used by sparsity-preservation assertions.  Every test
# that takes a sparse-aware path through the casadi NEP-13/18 bridge
# checks that the structural pattern is not silently densified.
_SP_LOWER3 = DM(Sparsity.lower(3), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])  # 6 nnz, 3x3
_SP_DIAG3  = DM(Sparsity.diag(3),  [1.0, 2.0, 3.0])                # 3 nnz, 3x3
_SP_LOWER4 = DM(Sparsity.lower(4),
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])  # 10 nnz, 4x4
_SP_VEC3   = DM(Sparsity.triplet(3, 1, [0, 2], [0, 0]),
                [1.5, 2.5])                                           # 2 nnz, 3x1


def _to_numpy(x):
    if isinstance(x, DM):
        return x.full()
    return np.asarray(x)


class _NumpyRefMixin(object):
    """Shared verification helpers."""

    def _check_nnz(self, label, np_func, *cas_args, expected_nnz, kwargs=None):
        """Call np_func through the casadi NEP-13/18 bridge on `cas_args`
        and assert the DM result has exactly `expected_nnz` structural
        nonzeros.  Catches handlers that quietly densify a sparse input
        (e.g. by going through np.asarray(x.full()))."""
        kwargs = kwargs or {}
        got = np_func(*cas_args, **kwargs)
        self.assertTrue(isinstance(got, DM),  # pyright: ignore[reportAttributeAccessIssue]
                        "[%s] got %s, expected DM" %
                        (label, type(got).__name__))
        self.assertEqual(got.nnz(), expected_nnz,  # pyright: ignore[reportAttributeAccessIssue]
                         "[%s] nnz=%d, expected %d (sparse input was densified?)" %
                         (label, got.nnz(), expected_nnz))

    def _close(self, got, ref, label, digits=11):
        """Compare casadi result `got` against numpy reference `ref`.

        The two array libraries disagree about 1-D shapes: numpy returns
        1-D from many reductions and elementwise ops on vectors, while
        casadi only has 2-D values.  We accept a numpy 1-D reference as
        equivalent to a column or row vector of the same length on the
        casadi side -- whichever orientation the casadi result has."""
        g = np.asarray(got.full() if isinstance(got, DM) else got)
        r = np.asarray(ref)
        if r.ndim == 0:
            r = r.reshape(1, 1)
        if g.ndim == 0:
            g = g.reshape(1, 1)
        if g.shape != r.shape:
            if r.ndim == 1 and g.ndim == 2 and g.size == r.size and \
                    (g.shape == (r.size, 1) or g.shape == (1, r.size)):
                r = r.reshape(g.shape)
        self.assertEqual(g.shape, r.shape,  # pyright: ignore[reportAttributeAccessIssue]
                         "[%s] shape: got %s, expected %s" %
                         (label, g.shape, r.shape))
        if g.size:
            d = np.abs(g.astype(float) - r.astype(float)).max()
            self.assertLess(d, 10 ** (-digits),  # pyright: ignore[reportAttributeAccessIssue]
                            "[%s] values: max abs diff %.3e\ngot=\n%s\nref=\n%s"
                            % (label, d, g, r))

    def _verify(self, label, np_func, *cas_args, digits=11,
                sym_types=(SX, MX), kwargs=None, np_args_override=None):
        """Run np_func through casadi dispatch on DM and each sym_type,
        and compare to the numpy reference computed on dense inputs.

        np_args_override: optional explicit args for the numpy reference
        call (used when numpy refuses 2-D vector inputs, e.g. np.cross
        of column vectors -- pass the flattened equivalents)."""
        kwargs = kwargs or {}

        # numpy reference
        if np_args_override is not None:
            np_args = np_args_override
        else:
            np_args = tuple(_to_numpy(a) if isinstance(a, DM) else a for a in cas_args)
        try:
            ref = np_func(*np_args, **kwargs)
        except Exception as e:
            self.fail("[%s] numpy reference itself raised: %s" % (label, e))  # pyright: ignore[reportAttributeAccessIssue]

        # DM
        try:
            got_dm = np_func(*cas_args, **kwargs)
        except Exception as e:
            self.fail("[%s] DM dispatch raised: %s" % (label, e))  # pyright: ignore[reportAttributeAccessIssue]
        self._close(got_dm, ref, label + " DM", digits=digits)

        # symbolic types
        for cls in sym_types:
            syms = []
            num_inputs = []
            sym_call_args = []
            for a in cas_args:
                if isinstance(a, DM):
                    s = cls.sym("x%d" % len(syms), a.shape[0], a.shape[1])
                    syms.append(s)
                    num_inputs.append(a)
                    sym_call_args.append(s)
                else:
                    sym_call_args.append(a)
            try:
                expr = np_func(*sym_call_args, **kwargs)
            except Exception as e:
                self.fail("[%s] %s dispatch raised: %s"  # pyright: ignore[reportAttributeAccessIssue]
                          % (label, cls.__name__, e))
            f = Function("f", syms, [expr])
            got = f(*num_inputs)
            self._close(got, ref, "%s %s" % (label, cls.__name__),
                        digits=digits)


class NumpyInteropTests(casadiTestCase, _NumpyRefMixin):
    """NEP-13 __array_ufunc__ and NEP-18 __array_function__ dispatch.

    Every value is compared to a pure-numpy reference; every test that
    feeds a structurally sparse DM also asserts the result's nnz, so a
    handler that quietly densifies (e.g. via np.asarray(x.full())) is
    caught regardless of whether the values still happen to match."""

    # -------------------- unary, all four shapes -------------------- #

    def _check_unary_all_shapes(self, np_func, x_factory=lambda x: x,
                                digits=11, sym_types=(SX, MX)):
        for label, x in _SHAPES_ELEMWISE:
            self._verify(np_func.__name__ + ":" + label,
                         np_func, x_factory(x),
                         digits=digits, sym_types=sym_types)

    def test_negative(self):    self._check_unary_all_shapes(np.negative)
    def test_positive(self):    self._check_unary_all_shapes(np.positive)
    def test_absolute(self):    self._check_unary_all_shapes(np.absolute, lambda x: -x)
    def test_fabs(self):        self._check_unary_all_shapes(np.fabs,     lambda x: -x)
    def test_square(self):      self._check_unary_all_shapes(np.square)
    def test_reciprocal(self):  self._check_unary_all_shapes(np.reciprocal)
    def test_sign(self):        self._check_unary_all_shapes(np.sign,     lambda x: x - 0.5)
    def test_sqrt(self):        self._check_unary_all_shapes(np.sqrt)
    def test_exp(self):         self._check_unary_all_shapes(np.exp)
    def test_expm1(self):       self._check_unary_all_shapes(np.expm1)
    def test_log(self):         self._check_unary_all_shapes(np.log)
    def test_log10(self):       self._check_unary_all_shapes(np.log10)
    def test_log1p(self):       self._check_unary_all_shapes(np.log1p)
    def test_sin(self):         self._check_unary_all_shapes(np.sin)
    def test_cos(self):         self._check_unary_all_shapes(np.cos)
    def test_tan(self):         self._check_unary_all_shapes(np.tan)
    def test_arcsin(self):      self._check_unary_all_shapes(np.arcsin)
    def test_arccos(self):      self._check_unary_all_shapes(np.arccos)
    def test_arctan(self):      self._check_unary_all_shapes(np.arctan)
    def test_sinh(self):        self._check_unary_all_shapes(np.sinh)
    def test_cosh(self):        self._check_unary_all_shapes(np.cosh)
    def test_tanh(self):        self._check_unary_all_shapes(np.tanh)
    def test_arcsinh(self):     self._check_unary_all_shapes(np.arcsinh)
    def test_arccosh(self):     self._check_unary_all_shapes(np.arccosh, lambda x: x + 1.0)
    def test_arctanh(self):     self._check_unary_all_shapes(np.arctanh)
    def test_floor(self):       self._check_unary_all_shapes(np.floor)
    def test_ceil(self):        self._check_unary_all_shapes(np.ceil)

    def test_logical_not(self):
        for label, x in [("vec", DM([0, 1, 0, 1])),
                         ("mat", DM([[0, 1], [1, 0]]))]:
            # casadi MX has no logic_not; SX is fine.
            self._verify("logical_not:" + label, np.logical_not, x,
                         sym_types=(SX,))

    # -------------------- binary, all four shapes ------------------- #

    def _check_binary_all_shapes(self, np_func, a_factory=lambda x: x,
                                 b_factory=lambda x: x + 0.1, digits=11,
                                 sym_types=(SX, MX)):
        for label, x in _SHAPES_ELEMWISE:
            a, b = a_factory(x), b_factory(x)
            self._verify(np_func.__name__ + ":" + label, np_func, a, b,
                         digits=digits, sym_types=sym_types)

    def test_add(self):         self._check_binary_all_shapes(np.add)
    def test_subtract(self):    self._check_binary_all_shapes(np.subtract)
    def test_multiply(self):    self._check_binary_all_shapes(np.multiply)
    def test_divide(self):      self._check_binary_all_shapes(np.divide)
    def test_true_divide(self): self._check_binary_all_shapes(np.true_divide)
    def test_power(self):       self._check_binary_all_shapes(np.power)
    def test_float_power(self): self._check_binary_all_shapes(np.float_power)
    def test_fmin(self):        self._check_binary_all_shapes(np.fmin)
    def test_fmax(self):        self._check_binary_all_shapes(np.fmax)
    def test_minimum(self):     self._check_binary_all_shapes(np.minimum)
    def test_maximum(self):     self._check_binary_all_shapes(np.maximum)
    def test_arctan2(self):     self._check_binary_all_shapes(np.arctan2)
    def test_hypot(self):       self._check_binary_all_shapes(np.hypot)
    def test_copysign(self):
        self._check_binary_all_shapes(np.copysign, b_factory=lambda x: x - 0.5)

    def test_compare(self):
        for op in (np.less, np.less_equal, np.greater, np.greater_equal,
                   np.equal, np.not_equal):
            for label, x in _SHAPES_ELEMWISE:
                self._verify(op.__name__ + ":" + label, op, x, x - 0.1)

    def test_logical_ops(self):
        for op in (np.logical_and, np.logical_or):
            self._verify(op.__name__,
                         op, DM([1, 0, 1, 0]), DM([1, 1, 0, 0]),
                         sym_types=(SX,))

    def test_mod_python_semantics(self):
        a = DM([-3.0, -7.0, -3.0, 3.0, 7.0, -3.5])
        b = DM([ 5.0,  5.0, -5.0, -5.0, -3.0, 2.0])
        self._verify("mod",       np.mod,       a, b)
        self._verify("remainder", np.remainder, a, b)

    def test_fmod_c_semantics(self):
        a = DM([-3.0, -7.0, 3.0, -3.0, 7.0, -3.5])
        b = DM([ 5.0,  5.0, -5.0, -5.0, -3.0, 2.0])
        self._verify("fmod", np.fmod, a, b)

    # -------------------- reductions ------------------------------- #

    def test_sum_global(self):
        for label, x in _SHAPES_ELEMWISE:
            self._verify("sum:" + label, np.sum, x)

    def test_sum_axis(self):
        # numpy axis=k drops that axis (shape collapses); casadi keeps
        # 2-D, so we use keepdims=True on the numpy side to align shapes.
        x = _MAT23
        for axis in (0, 1):
            self._verify("sum_axis%d" % axis, np.sum, x,
                         kwargs={"axis": axis, "keepdims": True})

    def test_cumsum_axis(self):
        # casadi.cumsum is per-axis; numpy default cumsum flattens.
        # We explicitly compare per-axis where they agree.
        for label, x in [("col", _COL), ("row", _ROW), ("mat", _MAT23)]:
            for axis in (0, 1):
                self._verify("cumsum:%s axis=%d" % (label, axis),
                             np.cumsum, x, kwargs={"axis": axis})

    def test_minmax_global(self):
        for label, x in _SHAPES_ELEMWISE:
            self._verify("max:" + label, np.max, x)
            self._verify("min:" + label, np.min, x)

    def test_all_any_global(self):
        # casadi MX has no logic_all/logic_any; SX is fine.
        self._verify("all_true",  np.all, DM([1, 1, 1]), sym_types=(SX,))
        self._verify("all_mixed", np.all, DM([1, 0, 1]), sym_types=(SX,))
        self._verify("any_false", np.any, DM([0, 0, 0]), sym_types=(SX,))
        self._verify("any_mixed", np.any, DM([0, 0, 1]), sym_types=(SX,))

    # -------------------- sparsity preservation -------------------- #

    def test_unary_preserves_sparsity(self):
        """Ufuncs with f(0)==0 must not densify a structurally sparse DM.

        A handler that internally goes through np.asarray(x.full()) would
        silently densify; this guards against that regression."""
        sp = Sparsity.lower(4)                # 10 structural nonzeros in 4x4
        x = DM(sp, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        x_nnz = x.nnz()
        zero_preserving = (np.negative, np.positive, np.square, np.sqrt,
                           np.sin, np.tan, np.arcsin, np.arctan,
                           np.sinh, np.tanh, np.arcsinh, np.arctanh,
                           np.expm1, np.log1p, np.floor, np.ceil,
                           np.fabs, np.absolute, np.sign)
        for f in zero_preserving:
            got = f(x)
            self.assertTrue(isinstance(got, DM),
                            "[%s] returned %s, not DM" % (f.__name__, type(got).__name__))
            self._check_nnz(f.__name__ + " sparse", f, x, expected_nnz=x_nnz)

    def test_binary_preserves_union_sparsity(self):
        """add/subtract on two sparse DMs must not produce a fully dense
        result when the union of the input patterns is itself sparse."""
        a = DM(Sparsity.lower(3), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])  # 6 nnz
        b = DM(Sparsity.diag(3),  [10.0, 20.0, 30.0])              # 3 nnz
        # union is just lower(3) since diag(3) is a subset; full would be 9.
        for f in (np.add, np.subtract):
            self._check_nnz(f.__name__ + " sparse", f, a, b,
                            expected_nnz=a.nnz())

    # ================== NEP-18 __array_function__ ================== #

    # -------------------- stacking --------------------------------- #

    def test_concatenate(self):
        a = DM([[1.0, 2.0], [3.0, 4.0]])
        b = DM([[5.0, 6.0], [7.0, 8.0]])
        for axis in (0, 1):
            self._verify("concat axis=%d" % axis,
                         lambda x, y, axis=axis: np.concatenate([x, y], axis=axis),
                         a, b)
        # sparse: concat preserves the sum of input nnz exactly
        for axis in (0, 1):
            self._check_nnz(
                "concat sparse axis=%d" % axis,
                lambda x, y, axis=axis: np.concatenate([x, y], axis=axis),
                _SP_LOWER3, _SP_DIAG3,
                expected_nnz=_SP_LOWER3.nnz() + _SP_DIAG3.nnz())

    def test_vstack(self):
        self._verify("vstack",
                     lambda x, y: np.vstack([x, y]),
                     DM([[1.0, 2.0], [3.0, 4.0]]),
                     DM([[5.0, 6.0]]))
        self._check_nnz("vstack sparse",
                        lambda x, y: np.vstack([x, y]),
                        _SP_LOWER3, _SP_DIAG3,
                        expected_nnz=_SP_LOWER3.nnz() + _SP_DIAG3.nnz())

    def test_hstack(self):
        self._verify("hstack",
                     lambda x, y: np.hstack([x, y]),
                     DM([[1.0, 2.0], [3.0, 4.0]]),
                     DM([[5.0], [6.0]]))
        self._check_nnz("hstack sparse",
                        lambda x, y: np.hstack([x, y]),
                        _SP_LOWER3, _SP_DIAG3,
                        expected_nnz=_SP_LOWER3.nnz() + _SP_DIAG3.nnz())

    # -------------------- shape ------------------------------------ #

    def test_reshape_C_order(self):
        v = DM([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        for shape in [(2, 3), (3, 2), (1, 6), (6, 1)]:
            self._verify("reshape_C %s" % (shape,),
                         lambda x, shape=shape: np.reshape(x, shape), v)
        # sparse: reshape changes the pattern but preserves nnz
        self._check_nnz("reshape_C sparse",
                        lambda x: np.reshape(x, (2, 8)),
                        _SP_LOWER4, expected_nnz=_SP_LOWER4.nnz())

    def test_reshape_F_order(self):
        v = DM([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        for shape in [(2, 3), (3, 2)]:
            self._verify("reshape_F %s" % (shape,),
                         lambda x, shape=shape: np.reshape(x, shape, order='F'),
                         v)
        self._check_nnz("reshape_F sparse",
                        lambda x: np.reshape(x, (2, 8), order='F'),
                        _SP_LOWER4, expected_nnz=_SP_LOWER4.nnz())

    def test_transpose(self):
        for label, x in _SHAPES_ELEMWISE:
            self._verify("transpose:" + label, np.transpose, x)
        # sparse: transpose preserves nnz exactly (transposes the pattern)
        self._check_nnz("transpose sparse", np.transpose,
                        _SP_LOWER4, expected_nnz=_SP_LOWER4.nnz())

    # -------------------- linear algebra --------------------------- #

    def test_dot_vec_vec(self):
        # numpy refuses to dot two (3,1) arrays -- compare against the
        # 1-D numpy result instead.
        v1, v2 = DM([1.0, 2.0, 3.0]), DM([4.0, 5.0, 6.0])
        self._verify("dot_vec_vec", np.dot, v1, v2,
                     np_args_override=(v1.full().ravel(), v2.full().ravel()))

    def test_dot_mat_mat(self):
        self._verify("dot_mat_mat", np.dot,
                     DM([[1.0, 2.0], [3.0, 4.0]]),
                     DM([[5.0, 6.0], [7.0, 8.0]]))

    def test_dot_mat_vec(self):
        self._verify("dot_mat_vec", np.dot,
                     DM([[1.0, 2.0], [3.0, 4.0]]),
                     DM([5.0, 6.0]))

    def test_matmul(self):
        self._verify("matmul", np.matmul,
                     DM([[1.0, 2.0], [3.0, 4.0]]),
                     DM([[5.0, 6.0], [7.0, 8.0]]))

    def test_inner_vec_vec(self):
        v1, v2 = DM([1.0, 2.0, 3.0]), DM([4.0, 5.0, 6.0])
        # numpy.inner of two (3,1) arrays returns (1,1) ((1,)*(1,) per row);
        # use ravelled refs for the conventional vector inner product.
        self._verify("inner_vec_vec", np.inner, v1, v2,
                     np_args_override=(v1.full().ravel(), v2.full().ravel()))

    def test_outer(self):
        self._verify("outer", np.outer,
                     DM([1.0, 2.0, 3.0]), DM([4.0, 5.0, 6.0]))
        # sparse: outer(a, b) has exactly nnz(a) * nnz(b) structural nonzeros
        sv2 = DM(Sparsity.triplet(3, 1, [1, 2], [0, 0]), [3.0, 4.0])
        self._check_nnz("outer sparse", np.outer, _SP_VEC3, sv2,
                        expected_nnz=_SP_VEC3.nnz() * sv2.nnz())

    def test_kron(self):
        self._verify("kron", np.kron,
                     DM([[1.0, 2.0]]), DM([[3.0], [4.0]]))
        # sparse: kron's nnz is the product of input nnz
        self._check_nnz("kron sparse", np.kron, _SP_LOWER3, _SP_DIAG3,
                        expected_nnz=_SP_LOWER3.nnz() * _SP_DIAG3.nnz())

    def test_cross(self):
        v1, v2 = DM([1.0, 2.0, 3.0]), DM([4.0, 5.0, 6.0])
        # numpy.cross requires last-axis length 2 or 3; (3,1) has last 1.
        # Compare against the ravelled-input numpy result.
        self._verify("cross", np.cross, v1, v2,
                     np_args_override=(v1.full().ravel(), v2.full().ravel()))

    def test_diag_vec_to_mat(self):
        # numpy.diag(1d) -> 2-D diag matrix.  A casadi column vector is
        # 2-D (n, 1); calling np.diag on it through casadi's bridge
        # returns the diagonal matrix.  Compare against numpy's 1-D
        # behavior directly.
        v = DM([1.0, 2.0, 3.0])
        self._verify("diag_vec_to_mat", np.diag, v,
                     np_args_override=(v.full().ravel(),))
        # sparse: diag of a sparse vector keeps only that vector's nnz
        # on the resulting diagonal
        self._check_nnz("diag_vec_to_mat sparse", np.diag, _SP_VEC3,
                        expected_nnz=_SP_VEC3.nnz())

    def test_diag_mat_to_vec(self):
        self._verify("diag_mat_to_vec", np.diag, _SPD3)
        # sparse: diag extracts the diagonal entries; lower(4) has a
        # fully populated diagonal -> 4 nnz
        self._check_nnz("diag_mat_to_vec sparse", np.diag, _SP_LOWER4,
                        expected_nnz=4)

    def test_trace(self):
        for label, x in [("2x2", DM([[1.0, 2.0], [3.0, 4.0]])),
                         ("spd3", _SPD3)]:
            self._verify("trace:" + label, np.trace, x)

    def test_linalg_norm_vector(self):
        v = DM([3.0, -4.0, 12.0])
        for ord in (None, 1, -1, 2, np.inf, -np.inf):
            self._verify("norm_vec ord=%s" % ord, np.linalg.norm, v,
                         kwargs={"ord": ord})

    def test_linalg_norm_matrix(self):
        A = DM([[1.0, 2.0, 3.0], [4.0, -5.0, 6.0]])
        for ord in (None, 'fro', 1, -1, 2, np.inf, -np.inf):
            # ord=2 is the spectral norm; no casadi equivalent for sym types.
            sym_types = () if ord == 2 else (SX, MX)
            self._verify("norm_mat ord=%s" % ord, np.linalg.norm, A,
                         kwargs={"ord": ord}, sym_types=sym_types)

    def test_linalg_solve_vec_rhs(self):
        self._verify("solve_vec", np.linalg.solve, _SPD3, _RHS3)

    def test_linalg_solve_mat_rhs(self):
        self._verify("solve_mat", np.linalg.solve, _SPD3, _RHS3M)

    def test_linalg_det(self):
        for label, A in [("2x2", DM([[1.0, 2.0], [3.0, 4.0]])),
                         ("spd3", _SPD3)]:
            # det(MX) is a known limitation in casadi core (no eval impl).
            self._verify("det:" + label, np.linalg.det, A,
                         sym_types=(SX,))

    def test_linalg_inv(self):
        for label, A in [("2x2", DM([[1.0, 2.0], [3.0, 4.0]])),
                         ("spd3", _SPD3)]:
            self._verify("inv:" + label, np.linalg.inv, A)

    def test_linalg_cholesky(self):
        # casadi.chol has no MX overload.
        self._verify("cholesky", np.linalg.cholesky, _SPD3,
                     sym_types=(SX,))

    # -------------------- repeat / tile ---------------------------- #

    def test_repeat_axis0(self):
        self._verify("repeat_axis0",
                     lambda x: np.repeat(x, 3, axis=0),
                     DM([[1.0, 2.0], [3.0, 4.0]]))
        # sparse: each row repeated 3 times -> nnz * 3
        self._check_nnz("repeat_axis0 sparse",
                        lambda x: np.repeat(x, 3, axis=0),
                        _SP_LOWER3, expected_nnz=_SP_LOWER3.nnz() * 3)

    def test_repeat_axis1(self):
        self._verify("repeat_axis1",
                     lambda x: np.repeat(x, 2, axis=1),
                     DM([[1.0, 2.0], [3.0, 4.0]]))
        self._check_nnz("repeat_axis1 sparse",
                        lambda x: np.repeat(x, 2, axis=1),
                        _SP_LOWER3, expected_nnz=_SP_LOWER3.nnz() * 2)

    def test_repeat_flat(self):
        for label, x in [("vec", DM([1.0, 2.0, 3.0])),
                         ("mat", DM([[1.0, 2.0], [3.0, 4.0]]))]:
            self._verify("repeat_flat:" + label,
                         lambda x: np.repeat(x, 2), x)
        # sparse: flattening repeat doubles the nnz
        self._check_nnz("repeat_flat sparse",
                        lambda x: np.repeat(x, 2),
                        _SP_LOWER3, expected_nnz=_SP_LOWER3.nnz() * 2)

    def test_tile_int(self):
        # Integer reps tiles along the last axis (numpy 2-D semantics).
        self._verify("tile_int", lambda x: np.tile(x, 2), DM([[1.0, 2.0]]))
        self._verify("tile_int_col", lambda x: np.tile(x, 2),
                     DM([[1.0], [2.0]]))
        # sparse: tile by k multiplies nnz by k
        self._check_nnz("tile_int sparse",
                        lambda x: np.tile(x, 2),
                        _SP_LOWER3, expected_nnz=_SP_LOWER3.nnz() * 2)

    def test_tile_2d(self):
        self._verify("tile_2d", lambda x: np.tile(x, (2, 3)),
                     DM([[1.0, 2.0]]))
        # sparse: tile by (m, n) multiplies nnz by m*n
        self._check_nnz("tile_2d sparse",
                        lambda x: np.tile(x, (2, 3)),
                        _SP_LOWER3, expected_nnz=_SP_LOWER3.nnz() * 6)

    # -------------------- like-constructors ------------------------ #

    def test_zeros_like(self):
        for label, x in _SHAPES_ELEMWISE:
            self._verify("zeros_like:" + label, np.zeros_like, x)

    def test_ones_like(self):
        for label, x in _SHAPES_ELEMWISE:
            self._verify("ones_like:" + label, np.ones_like, x)

    def test_full_like(self):
        for label, x in _SHAPES_ELEMWISE:
            self._verify("full_like:" + label,
                         lambda x: np.full_like(x, 7.5), x)

    # -------------------- conditional ------------------------------ #

    def test_where_elementwise(self):
        for label, x in [("col", _COL - 0.5), ("mat", _MAT23 - 0.5)]:
            self._verify("where:" + label,
                         lambda x: np.where(x > 0, x, -x), x)

    def test_where_three_arrays(self):
        cond = DM([1, 0, 1, 0])
        x    = DM([10.0, 20.0, 30.0, 40.0])
        y    = DM([-1.0, -2.0, -3.0, -4.0])
        self._verify("where_arrays", np.where, cond, x, y)

    # -------------------- range / clip / diff / roll --------------- #

    def test_linspace(self):
        # linspace endpoints are scalars; only test DM (numpy linspace
        # with int n always yields shape (n,)).
        for n in (2, 5, 10):
            ref = np.linspace(0.0, 1.0, n)
            got = np.linspace(DM(0), DM(1), n)
            self._close(got, ref, "linspace n=%d" % n, digits=12)

    def test_clip(self):
        self._verify("clip",
                     lambda x: np.clip(x, -1.0, 1.0),
                     DM([-3.0, -0.5, 0.0, 0.5, 3.0]))
        # sparse: bounds bracket 0 -> clip(0, -1, 1) = 0, sparsity preserved
        self._check_nnz("clip sparse",
                        lambda x: np.clip(x, -1.0, 1.0),
                        _SP_LOWER3, expected_nnz=_SP_LOWER3.nnz())

    def test_diff(self):
        col = DM([1.0, 3.0, 6.0, 10.0])
        row = DM([[1.0, 3.0, 6.0, 10.0]])
        mat = DM([[1.0, 2.0, 4.0, 7.0],
                  [0.0, 1.0, 3.0, 6.0],
                  [1.0, 1.0, 2.0, 4.0]])
        # default axis=-1: on a column the last axis is 1 (size 1),
        # giving an empty (n, 0) result -- matches numpy strictly.
        self._verify("diff_col_default", np.diff, col)
        # explicit axis along the populated dimension
        self._verify("diff_col_axis0", np.diff, col, kwargs={"axis": 0})
        self._verify("diff_row_axis1", np.diff, row, kwargs={"axis": 1})
        # higher-order on a row vector
        self._verify("diff_row_n2", np.diff, row, kwargs={"n": 2})
        self._verify("diff_row_n3", np.diff, row, kwargs={"n": 3})
        # matrix axis=0 and axis=1
        self._verify("diff_mat_axis0", np.diff, mat, kwargs={"axis": 0})
        self._verify("diff_mat_axis1", np.diff, mat, kwargs={"axis": 1})
        # n=0 is the no-op identity
        self._verify("diff_n0", np.diff, mat, kwargs={"n": 0})
        # repeated differencing on a matrix
        self._verify("diff_mat_n2_axis1", np.diff, mat,
                     kwargs={"n": 2, "axis": 1})

    def test_roll_1d(self):
        v = DM([1.0, 2.0, 3.0, 4.0, 5.0])
        for k in (-1, 0, 1, 2, 7):
            self._verify("roll_1d k=%d" % k,
                         lambda x: np.roll(x, k), v)
        # sparse: roll only shifts elements, never adds them
        for k in (-1, 0, 1, 2):
            self._check_nnz("roll_1d sparse k=%d" % k,
                            lambda x, k=k: np.roll(x, k),
                            _SP_VEC3, expected_nnz=_SP_VEC3.nnz())

    def test_roll_axis(self):
        m = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        for axis in (0, 1):
            for k in (-1, 0, 1, 4):
                self._verify("roll axis=%d k=%d" % (axis, k),
                             lambda x: np.roll(x, k, axis=axis), m)
        # sparse: same -- preserves nnz exactly along any axis
        for axis in (0, 1):
            for k in (-1, 0, 1):
                self._check_nnz(
                    "roll sparse axis=%d k=%d" % (axis, k),
                    lambda x, k=k, axis=axis: np.roll(x, k, axis=axis),
                    _SP_LOWER3, expected_nnz=_SP_LOWER3.nnz())

    # -------------------- atleast_*d -------------------------------- #

    def test_atleast(self):
        for label, x in _SHAPES_ELEMWISE:
            self._verify("atleast_1d:" + label, np.atleast_1d, x)
            self._verify("atleast_2d:" + label, np.atleast_2d, x)
        # sparse: atleast_*d on a 2-D sparse DM must be a structural
        # no-op (the matplotlib NEP-18 trap: pass-through must stay
        # sparse, not go through np.asarray(x.full())).
        for f in (np.atleast_1d, np.atleast_2d):
            self._check_nnz("%s sparse" % f.__name__, f,
                            _SP_LOWER4, expected_nnz=_SP_LOWER4.nnz())

    # -------------------- failure modes ---------------------------- #

    def test_unsupported_returns_clear_error(self):
        # Symbolic operations that have no casadi equivalent must
        # produce a clear failure rather than a misleading warning.
        with self.assertRaises(Exception):
            np.argmax(SX.sym("x", 3))


if __name__ == '__main__':
    unittest.main()
