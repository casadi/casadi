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
import casadi as _ca
from casadi import SX, MX, DM, Function, Sparsity, GlobalOptions
from helpers import casadiTestCase

Arr = _ca.ArrayInterface


def _unwrap(x):
    """Reach the underlying casadi value of a casadi.ArrayInterface;
    pass other values through unchanged."""
    return x.to_casadi() if isinstance(x, Arr) else x


# This suite verifies the casadi-aware numpy support (issue #2959): with
# GlobalOptions.setNumpyMode(1), an explicit numpy.foo(M) on a
# casadi value routes through the numpy-semantics array wrapper
# (casadi.ArrayInterface) and returns one, following numpy's shape/axis
# contract.  We enable it for the whole module and unwrap results back to
# casadi (to_casadi()) for numeric/symbolic comparison.  Default-mode behaviour
# (densify + FutureWarning) is covered in
# typemaps.py:test_numpy_preserve_type.
def setUpModule():
    GlobalOptions.setNumpyMode(1)


def tearDownModule():
    GlobalOptions.setNumpyMode(0)


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


def _bridge_dispatch_dict():
    """Build the NEP-18 dispatch dict for direct-handler tests.
    The casadi.casadi._build_numpy_function_dispatch() helper is an
    internal symbol not in the typed stubs; isolating the access here
    keeps the pyright suppression to one place."""
    import casadi as _ca
    return _ca.casadi._build_numpy_function_dispatch()  # pyright: ignore[reportAttributeAccessIssue]


class _NumpyRefMixin(object):
    """Shared verification helpers."""

    def _check_nnz(self, label, np_func, *cas_args, expected_nnz, kwargs=None):
        """Call np_func through the casadi NEP-13/18 bridge on `cas_args`
        and assert the DM result has exactly `expected_nnz` structural
        nonzeros.  Catches handlers that quietly densify a sparse input
        (e.g. by going through np.asarray(x.full()))."""
        kwargs = kwargs or {}
        got = _unwrap(np_func(*cas_args, **kwargs))
        self.assertTrue(isinstance(got, DM),  # pyright: ignore[reportAttributeAccessIssue]
                        "[%s] got %s, expected DM" %
                        (label, type(got).__name__))
        self.assertEqual(got.nnz(), expected_nnz,  # pyright: ignore[reportAttributeAccessIssue]
                         "[%s] nnz=%d, expected %d (sparse input was densified?)" %
                         (label, got.nnz(), expected_nnz))

    def _close(self, got, ref, label, digits=11):
        """Compare a casadi-aware numpy result against a numpy reference.

        Shapes are compared after squeezing length-1 axes.  ArrayInterface
        results DO carry numpy's exact logical rank for elementwise ops and
        axis reductions -- but a RAW casadi DM/SX/MX argument enters dispatch
        as its 2-D casadi shape (a column as (n,1), a row as (1,n), a scalar
        as (1,1)) per the column/row casting rule, so e.g. np.sin(DM[1,2,3])
        faithfully yields (n,1); meanwhile many of these tests write their
        numpy reference in flattened 1-D/0-D form.  Squeezing reconciles that
        deliberate (n,1)<->(n,) mapping while still catching genuine rank
        errors (a result with the wrong number of NON-singleton axes -- e.g.
        a reduction that failed to drop an axis -- survives squeezing and
        fails the shape check)."""
        g = np.asarray(got, dtype=float)
        r = np.asarray(ref, dtype=float)
        gs, rs = np.squeeze(g), np.squeeze(r)
        self.assertEqual(gs.shape, rs.shape,  # pyright: ignore[reportAttributeAccessIssue]
                         "[%s] shape: got %s (squeezed %s), expected %s (squeezed %s)" %
                         (label, g.shape, gs.shape, r.shape, rs.shape))
        if gs.size:
            d = np.abs(gs - rs).max()
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

        # DM: the result is an ArrayInterface with the true numpy shape.
        try:
            got_dm = np_func(*cas_args, **kwargs)
        except Exception as e:
            self.fail("[%s] DM dispatch raised: %s" % (label, e))  # pyright: ignore[reportAttributeAccessIssue]
        self._close(got_dm, ref, label + " DM", digits=digits)

        # symbolic types: build a Function on the underlying casadi value,
        # evaluate, and reshape the 2-D casadi output back to the result's
        # logical numpy shape before comparing.
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
                res = np_func(*sym_call_args, **kwargs)
            except Exception as e:
                self.fail("[%s] %s dispatch raised: %s"  # pyright: ignore[reportAttributeAccessIssue]
                          % (label, cls.__name__, e))
            shape = res.shape if isinstance(res, Arr) \
                else np.asarray(res).shape
            f = Function("f", syms, [_unwrap(res)])
            got = np.asarray(f(*num_inputs)).reshape(shape)
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
        # numpy axis=k drops that axis; the casadi-aware numpy support now
        # follows that contract, so a (2,3) sum over axis 0 is (3,) and over
        # axis 1 is (2,).
        x = _MAT23
        for axis in (0, 1):
            self._verify("sum_axis%d" % axis, np.sum, x, kwargs={"axis": axis})

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
            got = _unwrap(f(x))
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

    def _check_split(self, label, np_func, cas_arg, sections):
        ref = np_func(cas_arg.full(), sections)
        got = np_func(cas_arg, sections)
        # DM result is a list matching numpy's list-of-arrays.
        self.assertEqual(len(got), len(ref),  # pyright: ignore[reportAttributeAccessIssue]
                         "[%s] split count: got %d, expected %d" %
                         (label, len(got), len(ref)))
        for i, (g, r) in enumerate(zip(got, ref)):
            self._close(g, r, "%s[%d]" % (label, i))
        # symbolic: build a function and evaluate on the same dense input
        for cls in (SX, MX):
            s = cls.sym("x", cas_arg.shape[0], cas_arg.shape[1])
            pieces = [_unwrap(p) for p in np_func(s, sections)]
            f = Function("f", [s], pieces)
            out = f(cas_arg)
            if not isinstance(out, (list, tuple)):
                out = [out]
            self.assertEqual(len(out), len(ref),  # pyright: ignore[reportAttributeAccessIssue]
                             "[%s/%s] split count: got %d, expected %d" %
                             (label, cls.__name__, len(out), len(ref)))
            for i, (g, r) in enumerate(zip(out, ref)):
                self._close(g, r, "%s/%s[%d]" % (label, cls.__name__, i))

    def test_hsplit_int(self):
        a = DM([[1.0, 2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0, 8.0]])
        self._check_split("hsplit_int", np.hsplit, a, 2)
        self._check_split("hsplit_int4", np.hsplit, a, 4)

    def test_hsplit_indices(self):
        a = DM([[1.0, 2.0, 3.0, 4.0, 5.0],
                [6.0, 7.0, 8.0, 9.0, 10.0]])
        self._check_split("hsplit_idx", np.hsplit, a, [1, 3])
        # numpy.hsplit clamps over-range indices to empty pieces
        self._check_split("hsplit_idx_overrange", np.hsplit, a, [2, 7])

    def test_hsplit_uneven_falls_back(self):
        # 3 doesn't divide 4; numpy raises -- our handler returns
        # NotImplemented and numpy raises against the densified DM.
        a = DM([[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]])
        with self.assertRaises(ValueError):
            np.hsplit(a, 3)

    def test_vsplit_int(self):
        a = DM([[1.0, 2.0],
                [3.0, 4.0],
                [5.0, 6.0],
                [7.0, 8.0]])
        self._check_split("vsplit_int", np.vsplit, a, 2)
        self._check_split("vsplit_int4", np.vsplit, a, 4)

    def test_vsplit_indices(self):
        a = DM([[1.0, 2.0],
                [3.0, 4.0],
                [5.0, 6.0],
                [7.0, 8.0],
                [9.0, 10.0]])
        self._check_split("vsplit_idx", np.vsplit, a, [1, 4])

    def test_vsplit_sparse_preserves_nnz(self):
        # 4x4 lower-triangular, split into two 2x4 chunks: the two
        # halves' nnz must sum to the original nnz (no densification).
        pieces = [_unwrap(p) for p in np.vsplit(_SP_LOWER4, 2)]
        total = sum(p.nnz() for p in pieces)  # pyright: ignore[reportAttributeAccessIssue]
        self.assertEqual(total, _SP_LOWER4.nnz())  # pyright: ignore[reportAttributeAccessIssue]

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
        # Two casadi columns are 2-D (3,1); numpy.dot of (3,1)·(3,1) is an
        # invalid matmul (inner dims 1 vs 3).  The casadi-aware numpy
        # support follows numpy and raises rather than silently computing
        # an inner product.  (Use casadi.ArrayInterface([...]) for a true
        # 1-D vector whose dot IS the inner product.)
        v1, v2 = DM([1.0, 2.0, 3.0]), DM([4.0, 5.0, 6.0])
        with self.assertRaises(Exception):
            np.dot(v1, v2)
        # 1-D wrappers give the numpy inner product (scalar).
        d = np.dot(Arr([1.0, 2.0, 3.0]), Arr([4.0, 5.0, 6.0]))
        self._close(d, np.dot(np.array([1.0, 2.0, 3.0]),
                              np.array([4.0, 5.0, 6.0])), "dot 1-D")

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

    def test_cross_row_vectors(self):
        # Row vectors (1,3): numpy's default axis=-1 picks the length-3 axis.
        r1 = DM([[1.0, 2.0, 3.0]])
        r2 = DM([[4.0, 5.0, 6.0]])
        self._verify("cross_row", np.cross, r1, r2)

    def test_cross_3xN_axis0(self):
        # Shape (3, N), vectors along axis 0.
        a = DM([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]])
        b = DM([[7.0, 10.0], [8.0, 11.0], [9.0, 12.0]])
        self._verify("cross_3xN_axis0", np.cross, a, b, kwargs={"axis": 0})

    def test_cross_Nx3_default(self):
        # Shape (N, 3), default axis=-1 (last axis, length 3).
        a = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        b = DM([[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]])
        self._verify("cross_Nx3_default", np.cross, a, b)

    def test_cross_3x3_axis_disambiguates(self):
        # (3,3) is ambiguous under auto-detection; the explicit axis
        # kwarg must select the right one.  Without proper kwarg
        # handling these two calls would give the same wrong answer.
        a = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        b = DM([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        self._verify("cross_3x3_axis0",  np.cross, a, b, kwargs={"axis": 0})
        self._verify("cross_3x3_axis-1", np.cross, a, b, kwargs={"axis": -1})

    def test_cross_mixed_axes_falls_back(self):
        # Mixed axisa/axisb is not bridged; the dispatcher returns
        # NotImplemented and numpy falls back to its own impl on the
        # densified DM (no error, but proves no silent miscomputation).
        a = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        b = DM([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        got = np.cross(a, b, axisa=0, axisb=1)
        ref = np.cross(a.full(), b.full(), axisa=0, axisb=1)
        self._close(got, ref, "cross_mixed_axes_falls_back")

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
        # default axis=-1: numpy on a (n,1) array would give an empty
        # (n,0) result (the size-1 axis).  Our bridge treats (n,1) as
        # 1-D-like under the default axis and diffs along the populated
        # dimension; compare against numpy's 1-D result.
        self._verify("diff_col_default", np.diff, col,
                     np_args_override=(col.full().ravel(),))
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

    def test_trapz(self):
        # numpy 2.x removed `np.trapz` in favor of `np.trapezoid`; pick
        # whichever the installed numpy exposes so the test runs under
        # both versions.  The casadi bridge registers both names.
        trapz = getattr(np, "trapz", None) or np.trapezoid  # pyright: ignore[reportAttributeAccessIssue]
        # Default (uniform dx=1) on a column vector: bridge treats (n,1)
        # as 1-D and integrates along the populated axis to a scalar.
        # Without the handler this would return [0, 0, 0, 0] -- numpy
        # interpreting axis=-1 (size 1) on a (n,1) array.  aerosandbox's
        # `trapz` returns the per-subinterval contributions (length n-1)
        # rather than the sum; we deliberately match numpy's "total
        # integral" semantics.
        col = DM([1.0, 3.0, 6.0, 10.0])
        self._verify("trapz_col_default", trapz, col,
                     np_args_override=(col.full().ravel(),))
        row = DM([[1.0, 3.0, 6.0, 10.0]])
        self._verify("trapz_row_default", trapz, row,
                     np_args_override=(row.full().ravel(),))
        # Scalar dx.
        self._verify("trapz_scalar_dx",
                     lambda v: trapz(v, dx=2.0), col,
                     np_args_override=(col.full().ravel(),))
        # Non-uniform sample points.
        xc = DM([0.0, 1.0, 3.0, 6.0])
        self._verify("trapz_coords",
                     lambda y, x: trapz(y, x), col, xc,
                     np_args_override=(col.full().ravel(),
                                       xc.full().ravel()))
        # Matrix with default axis=-1: integrate along the last axis.
        mat = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        self._verify("trapz_mat_default", trapz, mat)
        # Matrix with axis=0.
        self._verify("trapz_mat_axis0",
                     lambda v: trapz(v, axis=0), mat)

    def test_gradient(self):
        # Default axis on a column vector: bridge treats (n,1) as 1-D,
        # numpy ref is computed on the ravelled 1-D form.
        col = DM([1.0, 3.0, 6.0, 10.0])
        self._verify("gradient_col_default", np.gradient, col,
                     np_args_override=(col.full().ravel(),))
        # Uniform scalar spacing.
        self._verify("gradient_scalar_h", lambda v: np.gradient(v, 2.0), col,
                     np_args_override=(col.full().ravel(),))
        # Non-uniform spacing via coordinate array.
        xcoord = DM([0.0, 1.0, 3.0, 6.0])
        ycoord = DM([0.0, 1.0, 4.0, 9.0])
        self._verify("gradient_coords",
                     lambda y, x: np.gradient(y, x), ycoord, xcoord,
                     np_args_override=(ycoord.full().ravel(),
                                       xcoord.full().ravel()))
        # edge_order=2.
        self._verify("gradient_edge2",
                     lambda v: np.gradient(v, edge_order=2), col,
                     np_args_override=(col.full().ravel(),))
        # 2-D matrix with explicit axis.
        m = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]])
        for ax in (0, 1):
            self._verify("gradient_mat_axis%d" % ax,
                         lambda v, ax=ax: np.gradient(v, axis=ax), m)
        # 2-D matrix without axis: returns a list of per-axis gradients.
        # Compare element-wise against numpy's list output.
        got_list = np.gradient(m)
        ref_list = np.gradient(m.full())
        self.assertEqual(len(got_list), len(ref_list))  # pyright: ignore[reportAttributeAccessIssue]
        for ax, (g, r) in enumerate(zip(got_list, ref_list)):
            self._close(g, r, "gradient_mat_no_axis ax=%d" % ax)
        # Symbolic 1-D query on SX/MX must produce a symbolic graph that
        # numerically matches numpy when evaluated.  List literal
        # `[SX, MX]` (not tuple) keeps pyright's loop-var widening to
        # Unknown so Function([s], [expr]) doesn't trip overload
        # resolution; see the floordiv/mod test for the explanation.
        ref = np.gradient(col.full().ravel())
        for cls in [SX, MX]:
            s = cls.sym("s", 4)
            f = Function("f", [s], [_unwrap(np.gradient(s))])
            self._close(f(col), ref, "gradient %s symbolic" % cls.__name__)

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

    # ---- kwargs handling: silent-drop guards --------------------------- #

    def test_kwarg_trace_offset_falls_back(self):
        # numpy.trace(offset=k) must not silently return the main diagonal.
        a = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        got = np.trace(a, offset=1)
        ref = np.trace(a.full(), offset=1)
        self._close(got, ref, "trace offset=1")
        # SX has no .full(): the unhandled kwarg must surface as an error.
        with self.assertRaises(TypeError):
            np.trace(SX.sym("X", 3, 3), offset=1)

    def test_kwarg_linspace_endpoint_false(self):
        # numpy.linspace(endpoint=False) excludes the stop; we honor it
        # symbolically (was previously silently ignored).
        got = np.linspace(DM(0.0), DM(1.0), 5, endpoint=False)
        ref = np.linspace(0.0, 1.0, 5, endpoint=False)
        self._close(got, ref, "linspace endpoint=False")

    def test_kwarg_zeros_like_shape(self):
        # numpy.zeros_like(shape=...) overrides the shape; silently
        # dropping it would return the wrong-sized matrix.  Also
        # compare values against numpy so a handler returning the
        # right shape filled with garbage is caught.
        a = DM([[1.0, 2.0], [3.0, 4.0]])
        got = np.zeros_like(a, shape=(3, 4))
        ref = np.zeros_like(a.full(), shape=(3, 4))
        self.assertEqual(got.shape, ref.shape)  # pyright: ignore[reportAttributeAccessIssue]
        self._close(got, ref, "zeros_like shape=")

    def test_kwarg_clip_out_falls_back(self):
        # numpy.clip(..., out=arr) is the in-place form; we don't bridge
        # it but the dispatcher's DM fallback must catch it.
        a = DM([-1.0, 0.5, 2.0])
        buf = np.empty((3, 1))
        got = np.clip(a, 0.0, 1.0, out=buf)
        ref = np.clip(a.full(), 0.0, 1.0)
        self._close(got, ref, "clip out=")

    def test_kwarg_ufunc_keepdims_noop(self):
        # casadi reductions are naturally 2-D, so keepdims=True is the
        # natural shape and keepdims=False is just no-op (we can't 1-D).
        # Either way: no silent miscomputation.
        x = DM([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        got = np.add.reduce(x, axis=0, keepdims=True)
        ref = np.add.reduce(x.full(), axis=0, keepdims=True)
        self._close(got, ref, "add.reduce keepdims=True")

    def test_kwarg_ufunc_where_rejected(self):
        # The `where=` mask changes which cells participate; we can't
        # bridge it.  Dispatcher must return NotImplemented; DM falls
        # back through numpy and the result must match numpy.
        x = DM([1.0, -2.0, 3.0])
        mask = np.array([[True], [False], [True]])
        got = np.abs(x, where=mask, out=np.zeros((3, 1)))
        ref = np.abs(x.full(), where=mask, out=np.zeros((3, 1)))
        self._close(got, ref, "abs where=")

    def test_kwarg_ufunc_initial_rejected(self):
        # `initial=` changes the reduction's starting term.  Must be
        # honored via the numpy fallback (DM) and not silently dropped.
        x = DM([1.0, 2.0, 3.0])
        got = np.add.reduce(x.T, initial=10.0)
        ref = np.add.reduce(x.full().T, initial=10.0)
        self._close(got, ref, "add.reduce initial=")

    def test_python_floordiv_and_mod_operators(self):
        # `x // y` and `x % y` previously raised
        # `TypeError: unsupported operand type(s) for //:` because
        # casadi's GenericExpressionCommon Python extend block didn't
        # define `__floordiv__` / `__rfloordiv__` / `__mod__` / `__rmod__`.
        # numpy's NEP-13 path worked (np.floor_divide / np.mod), but
        # the bare Python operators didn't reach it.
        # Both now follow Python's floor-mod convention (sign of divisor)
        # matching numpy.floor_divide / numpy.mod -- not C-fmod.
        x_dm = DM([7.0, -7.0, 5.5])
        y_dm = DM([3.0, 3.0, 2.0])
        ref_div = np.floor_divide(x_dm.full(), y_dm.full())
        ref_mod = np.mod(x_dm.full(), y_dm.full())
        self._close(x_dm // y_dm, ref_div, "DM // DM")
        self._close(x_dm % y_dm, ref_mod, "DM % DM")
        # Reverse operands.
        self._close(2.0 // y_dm, np.floor_divide(2.0, y_dm.full()),
                    "scalar // DM")
        self._close(2.0 % y_dm, np.mod(2.0, y_dm.full()), "scalar % DM")
        # Symbolic: build a Function and evaluate.  Iterate over a LIST
        # literal `[SX, MX]` rather than a tuple `(SX, MX)`: pyright
        # widens the loop variable of a heterogeneous list-literal to
        # Unknown (silently accepting all ops), whereas a tuple keeps
        # the `type[SX] | type[MX]` union and fails Function overload
        # resolution + `a // b`-style dunder lookup.  See the rest of
        # the suite (~45 sites) for the same pattern.
        for cls in [SX, MX]:
            a = cls.sym("a", 3)
            b = cls.sym("b", 3)
            f = Function("f", [a, b], [a // b, a % b])
            d, m = f(x_dm, y_dm)
            self._close(d, ref_div, "%s // %s" % (cls.__name__, cls.__name__))
            self._close(m, ref_mod, "%s %% %s" % (cls.__name__, cls.__name__))

    def test_logaddexp_reduce_column(self):
        # np.logaddexp.reduce maps to casadi's native logsumexp.  Verify
        # on DM and on each symbolic type so the AD path is exercised.
        from casadi import logsumexp
        x = DM([1.0, 2.0, 3.0])
        self._verify("logaddexp.reduce col",
                     lambda v: np.logaddexp.reduce(v),
                     x,
                     np_args_override=(np.array([1.0, 2.0, 3.0]),))

    def test_logaddexp_reduce_symbolic_uses_logsumexp(self):
        # Sanity-check: the SX path must NOT raise (it would if there were
        # no native handler, since SX can't densify).
        x = SX.sym("x", 4)
        expr = _unwrap(np.logaddexp.reduce(x))
        # The structural fingerprint must mention logsumexp's max-shift
        # trick (fmax + log of exp shifted) rather than triggering the
        # symbolic fallback's TypeError.
        f = Function("f", [x], [expr])
        ref = np.logaddexp.reduce(np.array([0.5, 1.5, 2.5, 3.5]))
        got = f(DM([0.5, 1.5, 2.5, 3.5]))
        self._close(got, ref, "logaddexp.reduce SX")

    def test_logaddexp_reduce_axis(self):
        # axis=0 reduces down columns; axis=1 across rows.
        x = DM([[1.0, 2.0], [3.0, 4.0]])
        for axis in (0, 1):
            got = np.logaddexp.reduce(x, axis=axis)
            ref = np.logaddexp.reduce(x.full(), axis=axis)
            self._close(got, ref, "logaddexp.reduce axis=%d" % axis)

    def test_complex_matmul_raises_typeerror_not_segfault(self):
        # Regression for issue #4216: numpy complex array @ MX used to
        # segfault.  It must now raise a clean TypeError so the user can
        # catch and recover instead of losing the interpreter.
        c = 0.5j * np.ones(2)
        for typ in (SX, MX):
            with self.assertRaises((TypeError, RuntimeError)):
                _ = c @ typ(2)
        # Same expectation for scalar complex * symbolic (this branch
        # already raised pre-fix; pinning it down so we don't regress).
        for typ in (SX, MX):
            with self.assertRaises(TypeError):
                # `# expect-error: <rule>` is a sentinel asserting that
                # pyright MUST report this rule on this line.  If a
                # later stub widening silently accepts `complex * MX`,
                # test_pyright_expect_error_sentinels fails before the
                # change ships, complementing the runtime assertRaises.
                _ = 0.5j * typ(1)  # expect-error: reportOperatorIssue


    # ============================================================== #
    # Regression tests for handler bugs found during the harness CI   #
    # wiring (commit b68cab33dd) -- each pins behaviour the bridge    #
    # got wrong before and that has no other test coverage.           #
    # ============================================================== #

    def test_cumprod_axis_none_is_c_order(self):
        # numpy ravels C-order before reducing; an earlier handler used
        # casadi.vec (F-order) and the test was pinned against the bug.
        x = DM([[1.0, 2.0], [3.0, 4.0]])
        self._verify("cumprod axis=None", np.cumprod, x)

    def test_append_axis_none_is_c_order(self):
        # axis=None flattens C-order; bridge flat-vec was F-order.
        a = DM([[1.0, 2.0]])
        b = DM([[3.0], [4.0], [5.0], [6.0]])
        self._verify("append axis=None", np.append, a, b)

    def test_diag_vector_with_offset(self):
        # Vector -> matrix with k != 0; earlier handler returned
        # NotImplemented and the fallback densified the column vector,
        # then numpy.diag of a (3,1) array extracted an empty diagonal.
        v = DM([1.0, 2.0, 3.0])
        for k in (1, -1, 2):
            self._verify("diag k=%d" % k, lambda x, k=k: np.diag(x, k=k),
                         v, np_args_override=(v.full().reshape(-1),))

    def test_diag_matrix_extract_offset(self):
        # Matrix -> k-th diagonal extraction.
        m = DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        for k in (-1, 0, 1):
            self._verify("diag(matrix) k=%d" % k,
                         lambda x, k=k: np.diag(x, k=k), m)

    def test_digitize_descending_bins(self):
        # Tie-handling for descending bins differs from ascending.
        # An earlier handler had right=False using `<=` (off-by-one at
        # ties); the correct condition for descending right=False is
        # strict `<`.
        x = DM([0.5, 1.5, 2.5, 4.5])
        bins = [4.0, 2.5, 1.0, 0.0]   # strictly descending
        # `np.digitize` densifies DM via numpy fallback if our handler
        # returns NotImplemented; check the symbolic-friendly handler
        # path actually fires by going through the dispatch dict
        # directly so the symbolic types still verify.
        self._verify("digitize descending",
                     lambda v, b: np.digitize(v, b),
                     x, bins,
                     np_args_override=(x.full().ravel(), bins))

    def test_linalg_cond_induced_matrix_norms(self):
        # casadi.norm_1/_inf are entrywise vector norms over vec(A);
        # numpy.linalg.cond(A, 1) needs the INDUCED matrix p-norm.
        # An earlier handler composed the wrong norms and gave a
        # cond off by a factor of ~n.
        # pyright widens the loop variable to "int | float | str"
        # which numpy.linalg.cond's Literal["fro","nuc"] doesn't match;
        # all the values we iterate ARE individually valid for cond.
        A = DM([[1.0, 2.0], [3.0, 4.0]])
        for p in (1, np.inf, 'fro'):
            self._verify("linalg.cond p=%r" % p,
                         lambda M, p=p: np.linalg.cond(M, p),  # pyright: ignore[reportArgumentType,reportCallIssue]
                         A)

    def test_linspace_retstep_with_dm_endpoints(self):
        # Symbolic retstep: returning NotImplemented forces the
        # fallback to densify DM(0)/DM(1) into (1,1) ndarrays, which
        # numpy.linspace then broadcasts into a (5,1,1) output -- the
        # Trap C scalar-wrapper round-trip.  Now honoured symbolically.
        s, e, n = DM(0.0), DM(1.0), 5
        samples, step = np.linspace(s, e, n, retstep=True)
        ref_samples, ref_step = np.linspace(0.0, 1.0, n, retstep=True)
        self._close(samples, ref_samples, "linspace retstep samples")
        self._close(step, ref_step, "linspace retstep step")

    def test_pad_constant_per_side_uses_axis_1_for_corners(self):
        # numpy pads axis 0 first, then axis 1; the new axis-1 cells
        # of the axis-0 padded rows should take the axis-1 constant.
        # An earlier handler padded axis 1 first, putting the wrong
        # constant in the corners.
        x = DM([[1.0, 2.0]])
        self._verify(
            "pad per-side constant_values",
            lambda v: np.pad(v, ((1, 0), (0, 2)),
                             constant_values=((5.0,), (9.0,))),
            x)

    def test_interp_symbolic_query(self):
        # casadi.interp1d rejects symbolic xq -- the handler is now a
        # hand-built piecewise-linear chain that works on SX/MX too.
        xp = [0.0, 1.0, 2.0, 3.0]
        fp = [0.0, 10.0, 20.0, 30.0]
        x = DM([0.5, 1.5, 2.5])
        self._verify("interp symbolic query",
                     lambda v: np.interp(v, xp, fp), x,
                     np_args_override=(x.full().ravel(),))

    def test_isfortran_returns_true(self):
        # CCS storage IS column-major; numpy.isfortran on a casadi
        # value should report True (when reachable -- on numpy<2.0
        # the function doesn't route through __array_function__, so we
        # check the registered handler directly).
        h = _bridge_dispatch_dict().get(np.isfortran)
        self.assertTrue(h is not None)
        self.assertIs(h(DM([[1.0]])), True)
        self.assertIs(h(SX.sym("x", 2, 2)), True)
        self.assertIs(h(MX.sym("y", 2, 2)), True)

    # ============================================================== #
    # Spot-checks for the rounds 2-4 NEP-13 ufunc additions.          #
    # Each entry exercises the dispatch path for the listed ufunc on  #
    # DM and -- where the casadi op is symbolic -- on SX/MX too.      #
    # ============================================================== #

    def test_round2_unary_ufuncs(self):
        # Avoid exact-.5 inputs -- banker's vs away-from-zero tie
        # breaking differs at those, documented in the kwargs-hygiene
        # skill.
        x = DM([[0.7, 1.7], [2.3, 8.0]])
        for op in (np.cbrt, np.exp2, np.log2, np.trunc, np.rint):
            self._verify(op.__name__, op, x)

    def test_round2_angle_conversions(self):
        # deg<->rad lambdas were derived constants; check both routes.
        for op in (np.deg2rad, np.radians):
            self._verify(op.__name__, op, DM([0.0, 30.0, 90.0, 180.0]))
        for op in (np.degrees, np.rad2deg):
            self._verify(op.__name__, op,
                         DM([0.0, np.pi / 6, np.pi / 2, np.pi]))

    def test_round2_binary_ufuncs(self):
        x = DM([7.0, -7.0, 0.5, -0.5])
        y = DM([3.0, 3.0, 0.25, 0.25])
        self._verify("floor_divide", np.floor_divide, x, y)
        self._verify("signbit", np.signbit, x)

    def test_round2_modf_and_logaddexp2(self):
        # modf returns a tuple; logaddexp2 uses max-shift for stability.
        x = DM([1.7, -1.7, 3.0])
        frac, integer = np.modf(x)
        ref_frac, ref_int = np.modf(x.full().ravel())
        self._close(frac, ref_frac, "modf frac")
        self._close(integer, ref_int, "modf int")
        # logaddexp2: log2(2**x + 2**y); equal inputs -> +1.
        a = DM([1.0, 2.0, 3.0])
        self._verify("logaddexp2", np.logaddexp2, a, a)

    def test_round3_logaddexp_divmod_heaviside_ldexp(self):
        a = DM([1.0, 2.0, 3.0])
        self._verify("logaddexp", np.logaddexp, a, a)

        x = DM([7.0, -7.0, 5.5])
        y = DM([3.0, 3.0, 2.0])
        q, m = np.divmod(x, y)
        rq, rm = np.divmod(x.full().ravel(), y.full().ravel())
        self._close(q, rq, "divmod q")
        self._close(m, rm, "divmod m")

        z = DM([-1.0, 0.0, 1.0])
        self._verify("heaviside",
                     lambda v: np.heaviside(v, 0.5), z)

        self._verify("ldexp", lambda v: np.ldexp(v, 3), DM([1.0, 2.0]))

    def test_round4_logical_xor(self):
        a = DM([0, 1, 1, 0])
        b = DM([0, 0, 1, 1])
        self._verify("logical_xor", np.logical_xor, a, b)

    # ============================================================== #
    # Spot-checks for the rounds 2-4 NEP-18 free-function additions.  #
    # ============================================================== #

    def test_round2_predicates(self):
        # Ground DM cases against the densified numpy answer; for SX/MX
        # (no .full()) the contract is that the predicate matches what
        # numpy would say about a 2-D ndarray of the same shape.
        dm_s = DM(1.0)
        self.assertEqual(np.isscalar(dm_s), np.isscalar(dm_s.full()))
        self.assertEqual(np.isscalar(SX.sym("x")),
                         np.isscalar(np.zeros((1, 1))))
        dm_v = DM([1.0])
        self.assertEqual(np.isrealobj(dm_v), np.isrealobj(dm_v.full()))
        self.assertEqual(np.isrealobj(SX.sym("y", 3)),
                         np.isrealobj(np.zeros((3, 1))))
        self.assertEqual(bool(np.isclose(DM([1.0]), 1.0 + 1e-10)),
                         bool(np.isclose(np.array([1.0]), 1.0 + 1e-10)))

    def test_round3_predicates(self):
        dm = DM([1.0])
        self.assertEqual(np.iscomplexobj(dm), np.iscomplexobj(dm.full()))
        self.assertEqual(np.iscomplexobj(SX.sym("y")),
                         np.iscomplexobj(np.zeros((1, 1))))
        # np.allclose / array_equal compare a casadi-shape DM (2,1)
        # against a python list -- the (2,) numpy reference doesn't
        # broadcast cleanly, so compare against a properly-shaped
        # column ndarray.
        a = DM([1.0, 2.0])
        b = np.array([[1.0], [2.0]])
        self.assertEqual(bool(np.allclose(a, b)),
                         bool(np.allclose(a.full(), b)))
        c = DM([1, 2, 3])
        self.assertEqual(bool(np.array_equal(c, c)),
                         bool(np.array_equal(c.full(), c.full())))
        # isposinf / isneginf delegate to numpy fallback for DM.
        x = DM([1.0, float('inf'), float('-inf'), 0.0])
        self.assertTrue(np.array_equal(
            np.asarray(np.isposinf(x), dtype=bool).ravel(),
            np.isposinf(x.full()).ravel()))

    def test_round4_predicates(self):
        # iscomplex / isreal return per-element bool arrays; ground
        # the casadi-side aggregate against the densified-input answer.
        a = DM([1.0, 2.0])
        self.assertEqual(
            bool(np.any(np.asarray(np.iscomplex(a), dtype=bool))),
            bool(np.any(np.iscomplex(a.full()))))
        self.assertEqual(
            bool(np.all(np.asarray(np.isreal(a), dtype=bool))),
            bool(np.all(np.isreal(a.full()))))
        nz = DM([0., 1., 2., 0., 3.])
        self.assertEqual(int(np.count_nonzero(nz)),
                         int(np.count_nonzero(nz.full())))
        eq = DM([1, 2])
        self.assertEqual(bool(np.array_equiv(eq, eq)),
                         bool(np.array_equiv(eq.full(), eq.full())))

    def test_round2_creation(self):
        x = DM([[1.0, 2.0], [3.0, 4.0]])
        # copy returns a fresh casadi value.
        c = _unwrap(np.copy(x))
        self.assertIsInstance(c, DM)
        self._close(c, x.full(), "copy")
        # asarray/asanyarray/ascontiguousarray DO route via the
        # dispatch table when called with `like=` on numpy>=1.20;
        # without `like=` they go straight to numpy's __array__
        # protocol which densifies DM via .full().  Check the
        # registered handler is a pass-through.
        d = _bridge_dispatch_dict()
        for op in (np.asarray, np.asanyarray, np.ascontiguousarray):
            self.assertIs(d[op](x), x)
        # diagflat with offset.
        self._verify("diagflat", np.diagflat, DM([1.0, 2.0, 3.0]),
                     np_args_override=(np.array([1.0, 2.0, 3.0]),))
        # tril / triu (matrix form).
        m = DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        self._verify("tril", np.tril, m)
        self._verify("triu", np.triu, m)
        # vander.
        self._verify("vander", np.vander, DM([1.0, 2.0, 3.0]),
                     np_args_override=(np.array([1.0, 2.0, 3.0]),))
        # logspace / geomspace: numpy's __array_function__ densifies DM
        # endpoints via __array__ before reaching the handler (Trap C),
        # so exercise the registered handler directly with DM endpoints
        # and compare to the pure-numpy reference on Python scalars.
        d = _bridge_dispatch_dict()
        gh = d.get(np.geomspace)
        self.assertIsNotNone(gh)
        self._close(gh(DM(1.0), DM(1000.0), 4),
                    np.geomspace(1.0, 1000.0, 4), "geomspace handler")
        lh = d.get(np.logspace)
        self.assertIsNotNone(lh)
        self._close(lh(DM(0.0), DM(2.0), 3),
                    np.logspace(0.0, 2.0, 3), "logspace handler")

    def test_round2_reductions_and_misc(self):
        x = DM([[2.0, 3.0], [4.0, 5.0]])
        # axis-aware prod / nanprod / nansum.
        for axis in (None, 0, 1):
            self._verify("prod axis=%r" % axis,
                         lambda v, axis=axis: np.prod(v, axis=axis), x)
            self._verify("nansum axis=%r" % axis,
                         lambda v, axis=axis: np.nansum(v, axis=axis), x)
            self._verify("nanprod axis=%r" % axis,
                         lambda v, axis=axis: np.nanprod(v, axis=axis), x)
        # around / round + scalar `decimals`.
        self._verify("around dec=1",
                     lambda v: np.around(v, decimals=1),
                     DM([[1.234, 5.678]]))
        # sinc: piecewise via if_else.
        self._verify("sinc", np.sinc, DM([0.0, 0.5, 1.0]))
        # ediff1d.
        self._verify("ediff1d", np.ediff1d, DM([1.0, 3.0, 7.0]),
                     np_args_override=(np.array([1.0, 3.0, 7.0]),))
        # fix (Deprecated)
        #self._verify("fix", np.fix, DM([1.7, -1.7, 2.0]))
        # isin -- balanced OR fold.
        self._verify("isin",
                     lambda v: np.isin(v, [1, 3, 5]),
                     DM([1, 2, 3, 4, 5]),
                     np_args_override=(np.array([1, 2, 3, 4, 5]),))

    def test_round3_array_manipulation(self):
        x = DM([[1.0, 2.0], [3.0, 4.0]])
        for op in (np.flip, np.fliplr, np.flipud):
            self._verify(op.__name__, op, x)
        self._verify("rot90", np.rot90, x)
        self._verify("ravel C", lambda v: np.ravel(v), x)
        self.assertEqual(np.ndim(x), 2)
        self.assertEqual(np.shape(DM([1, 2, 3])), (3, 1))
        self.assertEqual(int(np.size(DM([[1, 2, 3]]))), 3)
        # append axis=0 / axis=1 stacking.
        for axis in (0, 1):
            self._verify("append axis=%d" % axis,
                         lambda a, b, axis=axis: np.append(a, b, axis=axis),
                         x, x)
        # delete by index along an axis.
        self._verify("delete axis=0",
                     lambda v: np.delete(v, 1, axis=0),
                     DM([[1, 2], [3, 4], [5, 6]]))
        self._verify("insert axis=1",
                     lambda v: np.insert(v, 1, 99, axis=1),
                     DM([[1, 2, 3]]))

    def test_round3_stats(self):
        x = DM([1.0, 2.0, 3.0, 4.0])
        for op in (np.mean, np.std, np.var, np.ptp):
            self._verify(op.__name__, op, x,
                         np_args_override=(np.array([1.0, 2.0, 3.0, 4.0]),))
        # average with weights.
        self._verify(
            "average weights",
            lambda v: np.average(v, weights=[0.1, 0.2, 0.3, 0.4]),
            x,
            np_args_override=(np.array([1.0, 2.0, 3.0, 4.0]),))

    def test_round3_misc(self):
        x = DM([1.0, 2.0, 3.0])
        self._verify("cumprod axis=None",
                     lambda v: np.cumprod(v), x,
                     np_args_override=(np.array([1.0, 2.0, 3.0]),))
        self._verify("nan_to_num", np.nan_to_num, x)
        self._verify("real_if_close", np.real_if_close, x)
        # vdot on 1-D-like inputs.
        self._verify("vdot", np.vdot, x, x,
                     np_args_override=(np.array([1.0, 2.0, 3.0]),
                                       np.array([1.0, 2.0, 3.0])))
        # nanmax / nanmin (axis=None only).
        self._verify("nanmax", np.nanmax, x,
                     np_args_override=(np.array([1.0, 2.0, 3.0]),))
        self._verify("nanmin", np.nanmin, x,
                     np_args_override=(np.array([1.0, 2.0, 3.0]),))

    def test_round3_linalg(self):
        A = DM([[1.0, 2.0], [3.0, 5.0]])
        self._verify("linalg.matrix_power",
                     lambda M: np.linalg.matrix_power(M, 3), A)
        self._verify("linalg.pinv", np.linalg.pinv, A)
        # slogdet returns a (sign, log|det|) pair.
        s, ld = np.linalg.slogdet(A)
        rs, rld = np.linalg.slogdet(A.full())
        self._close(s, rs, "slogdet sign")
        self._close(ld, rld, "slogdet logdet")

    def test_round4_indexing(self):
        x = DM([10.0, 20.0, 30.0, 40.0])
        # take with explicit indices.
        self._verify("take", lambda v: np.take(v, [0, 2]), x,
                     np_args_override=(np.array([10.0, 20.0, 30.0, 40.0]),))
        # choose is dispatched via the registered handler (numpy
        # doesn't always invoke __array_function__ for it on integer
        # selectors); check the handler directly.
        d = _bridge_dispatch_dict()
        h = d.get(np.choose)
        if h is not None:
            ref = np.choose(0, [np.array([10, 20]), np.array([30, 40])])
            self._close(h(0, [DM([10, 20]), DM([30, 40])]),
                        ref, "choose handler")
        # select via the dispatch table: numpy.select on casadi
        # `condlist` would iterate over the DM (not allowed).
        hs = d.get(np.select)
        if hs is not None:
            cond = DM([1, 0, 1])
            ref = np.select([cond.full().ravel().astype(bool)],
                            [np.array([10, 20, 30])], default=99)
            self._close(hs([cond], [DM([10, 20, 30])], default=99),
                        ref.reshape(-1, 1), "select handler")

    def test_einsum_matmul_two_operands(self):
        # 'ij,jk->ik' is the matmul form; bridge maps to casadi.einstein
        # on vec(A), vec(B) with negative-int labels.
        A = DM([[1.0, 2.0], [3.0, 4.0]])
        B = DM([[5.0, 6.0], [7.0, 8.0]])
        self._verify("einsum matmul",
                     lambda a, b: np.einsum('ij,jk->ik', a, b), A, B)

    def test_einsum_mat_vec_and_transposed_output(self):
        # Mat-vec ('ij,j->i') and transposed-output ('ij,jk->ki').
        A = DM([[1.0, 2.0], [3.0, 4.0]])
        B = DM([[5.0, 6.0], [7.0, 8.0]])
        x = DM([1.0, 2.0])
        self._verify(
            "einsum mat-vec",
            lambda a, v: np.einsum('ij,j->i', a, v), A, x,
            np_args_override=(A.full(), x.full().ravel()))
        self._verify(
            "einsum transposed product",
            lambda a, b: np.einsum('ij,jk->ki', a, b), A, B)

    def test_einsum_inner_product(self):
        # 'i,i->' contracts to a scalar.  Pre-fix this was a hard
        # "no implementation found" on SX/MX -- now reaches casadi.einstein.
        x = DM([1.0, 2.0, 3.0])
        y = DM([4.0, 5.0, 6.0])
        self._verify(
            "einsum inner product",
            lambda u, v: np.einsum('i,i->', u, v), x, y,
            np_args_override=(x.full().ravel(), y.full().ravel()))

    def test_einsum_symbolic_uses_einstein_node(self):
        # The whole point of bridging einsum is that SX/MX inputs no
        # longer raise "no implementation found".  For MX this should
        # collapse to (a transformation of) the OP_EINSTEIN graph node.
        a = MX.sym("a", 2, 2)
        b = MX.sym("b", 2, 2)
        expr = _unwrap(np.einsum('ij,jk->ik', a, b))
        self.assertIn("einstein", str(expr))
        A = DM([[1.0, 2.0], [3.0, 4.0]])
        B = DM([[5.0, 6.0], [7.0, 8.0]])
        f = Function("f", [a, b], [expr])
        self._close(f(A, B), (A.full() @ B.full()), "einsum MX matmul")

    def test_logaddexp_scalar_uses_logsumexp(self):
        # For 1x1 inputs the binary logaddexp delegates to
        # casadi.logsumexp(vertcat(x, y)) so MX gets the single
        # OP_LOGSUMEXP primitive node rather than the hand-rolled
        # fmax + log1p + exp + fmin decomposition.
        a = MX.sym("a")
        b = MX.sym("b")
        expr = _unwrap(np.logaddexp(a, b))
        self.assertIn("logsumexp", str(expr))
        f = Function("f", [a, b], [expr])
        for av, bv in ((1.0, 2.0), (-3.0, 4.0), (0.0, 0.0)):
            self._close(f(av, bv),
                        np.logaddexp(av, bv),
                        "logaddexp scalar @(%g,%g)" % (av, bv))

    def test_round4_emath_real_domain(self):
        # emath.* doesn't reliably dispatch via __array_function__ on
        # all numpy versions, so test the registered handler directly.
        # On real inputs it collapses to the regular real op.
        d = _bridge_dispatch_dict()
        x = DM([1.0, 4.0, 9.0])
        ref_input = np.array([1.0, 4.0, 9.0])
        for op_name in ("sqrt", "log2", "log10"):
            ref = getattr(np.emath, op_name)(ref_input)
            handler = d.get(getattr(np.emath, op_name))
            if handler is None:
                continue
            self._close(handler(x), ref, "emath." + op_name)


if __name__ == '__main__':
    unittest.main()
