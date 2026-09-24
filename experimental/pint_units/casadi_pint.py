"""Proof of concept: casadi expressions with pint units.

`CQ` bundles a casadi expression (SX/MX/DM) with a pint unit by holding a
`pint.Quantity` whose magnitude is the casadi expression.  All unit algebra
(multiplication, conversion, dimension checks) is delegated to pint; `CQ`
only adds what pint lacks for casadi:

  * `__SX__` / `__MX__` / `__DM__`: casadi's conversion hooks, so a `CQ`
    can be passed straight to casadi (ca.Function, ca.vertcat, nlpsol, ...).
    The expression is handed over in the CQ's own unit.  With
    `CQ.strict = True` only dimensionless quantities may cross into casadi.
  * No `__getattr__` forwarding.  A raw `pint.Quantity` forwards SWIG's
    `.this` pointer to its magnitude, so `ca.sin(3*ureg.m)` silently strips
    the unit; `CQ` does not.
  * casadi-flavoured helpers that keep track of units: `vertcat`,
    `jacobian`, `mtimes`, `sin`/`exp`/... (dimensionless only),
    `sqrt`, `sumsqr`.

    ureg = pint.UnitRegistry()
    x = CQ.sym("x", ureg.km)         # SX symbol, unit km
    v = CQ.sym("v", "m/s", kind=ca.MX)
    t = x / v                        # CQ in km*s/m
    t.to("s")                        # CQ, expression 1000*x/v
"""

import casadi as ca
import pint

_CASADI_TYPES = (ca.SX, ca.MX, ca.DM)


class CQ(object):
    """casadi expression + pint unit."""

    # Win mixed binary ops against numpy arrays and casadi types
    # (casadi's operators yield to a higher __array_priority__).
    __array_priority__ = 2000.0

    # If True, __SX__/__MX__/__DM__ refuse dimensional quantities.
    strict = False

    ureg = None  # default registry, set with CQ.set_registry(...)

    def __init__(self, expr, unit=None, ureg=None):
        if isinstance(expr, CQ):
            q = expr.q if unit is None else expr.q.to(unit)
        elif isinstance(expr, pint.Quantity):
            q = expr if unit is None else expr.to(unit)
        else:
            ureg = ureg or CQ.ureg
            if ureg is None:
                raise ValueError("no unit registry: pass ureg= or call CQ.set_registry")
            q = ureg.Quantity(expr, "dimensionless" if unit is None else unit)
        self.q = q

    @classmethod
    def set_registry(cls, ureg):
        cls.ureg = ureg

    @classmethod
    def sym(cls, name, unit, *shape, kind=ca.SX, ureg=None):
        return cls(kind.sym(name, *shape), unit, ureg=ureg)

    # -- basic accessors ---------------------------------------------------
    @property
    def expr(self):
        """The casadi expression, in `self.units`."""
        return self.q.magnitude

    m = magnitude = expr

    @property
    def units(self):
        return self.q.units

    @property
    def dimensionless(self):
        return self.q.dimensionless

    @property
    def shape(self):
        return self.expr.shape

    @property
    def T(self):
        return CQ(self.q.magnitude.T * self.q.units)

    def to(self, unit):
        return CQ(self.q.to(unit))

    def to_base_units(self):
        return CQ(self.q.to_base_units())

    def m_as(self, unit):
        """casadi expression expressed in `unit` (explicit unit stripping)."""
        return self.q.to(unit).magnitude

    def __repr__(self):
        return "CQ(%s, '%s')" % (self.expr, self.units)

    __str__ = __repr__

    # -- casadi conversion hooks -------------------------------------------
    def _to_casadi(self, kind):
        q = self.q
        if not q.dimensionless:
            if CQ.strict:
                raise pint.DimensionalityError(
                    q.units, "dimensionless",
                    extra_msg=" (CQ.strict: convert explicitly with .m_as(unit))")
        else:
            q = q.to("dimensionless")   # e.g. m/km -> factor 1e-3 folded in
        e = q.magnitude
        if kind is ca.DM and not isinstance(e, ca.DM):
            if isinstance(e, (ca.SX, ca.MX)):
                if not e.is_constant():
                    return None         # symbolic: not a DM
                e = ca.evalf(e)
        return kind(e)

    def __SX__(self):
        return None if isinstance(self.expr, ca.MX) else self._to_casadi(ca.SX)

    def __MX__(self):
        return None if isinstance(self.expr, ca.SX) else self._to_casadi(ca.MX)

    def __DM__(self):
        return self._to_casadi(ca.DM)

    def __float__(self):
        return float(self.__DM__())

    # -- arithmetic: delegate to pint --------------------------------------
    @staticmethod
    def _q(o):
        if isinstance(o, CQ):
            return o.q
        if isinstance(o, (pint.Unit,)):
            return 1 * o
        return o

    @staticmethod
    def _wrap(r):
        return CQ(r) if isinstance(r, pint.Quantity) else r

    def _binop(op):
        def f(self, other):
            return CQ._wrap(getattr(self.q, op)(CQ._q(other)))
        f.__name__ = op
        return f

    for _op in ("add", "sub", "mul", "truediv", "pow"):
        locals()["__%s__" % _op] = _binop("__%s__" % _op)
        locals()["__r%s__" % _op] = _binop("__r%s__" % _op)
    del _op, _binop

    # matrix product: pint would route through np.matmul, which cannot mix
    # a numpy magnitude with a casadi one -> multiply magnitudes with casadi
    @staticmethod
    def _matmul(a, b):
        a, b = CQ._q(a), CQ._q(b)
        ma, ua = (a.magnitude, a.units) if isinstance(a, pint.Quantity) else (a, 1)
        mb, ub = (b.magnitude, b.units) if isinstance(b, pint.Quantity) else (b, 1)
        return CQ(ca.mtimes(ma, mb) * (ua * ub))

    def __matmul__(self, other):
        return CQ._matmul(self, other)

    def __rmatmul__(self, other):
        return CQ._matmul(other, self)

    # comparisons return plain casadi expressions (constraints) with both
    # sides converted to self.units.  (pint's own __eq__ calls bool() on the
    # magnitude, which symbolic casadi expressions refuse.)
    def _cmp(op):
        def f(self, other):
            o = CQ._q(other)
            if isinstance(o, pint.Quantity):
                o = o.to(self.units).magnitude
            elif not self.dimensionless:
                raise pint.DimensionalityError(self.units, "dimensionless")
            return getattr(self.expr, op)(o)
        f.__name__ = op
        return f

    for _op in ("lt", "le", "gt", "ge", "eq", "ne"):
        locals()["__%s__" % _op] = _cmp("__%s__" % _op)
    del _op, _cmp
    __hash__ = object.__hash__

    def __neg__(self):
        return CQ(-self.q)

    def __pos__(self):
        return self

    def __abs__(self):
        return CQ(fabs(self))

    def __getitem__(self, k):
        return CQ(self.expr[k] * self.units)

    # numpy ufuncs (np.sin(cq), ...): let pint do the unit bookkeeping
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if ufunc.__name__ == "matmul" and method == "__call__":
            return CQ._matmul(*inputs)
        r = getattr(ufunc, method)(*[CQ._q(i) for i in inputs], **kwargs)
        return CQ._wrap(r)


# -- unit-aware casadi helpers ---------------------------------------------

def _cq(x):
    return x if isinstance(x, CQ) else CQ(x)


def _dimless(x, fname):
    x = _cq(x)
    if not x.dimensionless:
        raise pint.DimensionalityError(x.units, "dimensionless",
                                       extra_msg=" in %s()" % fname)
    return x.m_as("dimensionless")


def _unary_dimless(fname):
    cafun = getattr(ca, fname)

    def f(x):
        return CQ(cafun(_dimless(x, fname)), "dimensionless", ureg=_cq(x).q._REGISTRY)
    f.__name__ = fname
    return f


sin, cos, tan, exp, log, tanh, atan, asin, acos = [
    _unary_dimless(n) for n in ("sin", "cos", "tan", "exp", "log", "tanh",
                                "atan", "asin", "acos")]


def sqrt(x):
    x = _cq(x)
    return CQ(ca.sqrt(x.expr) * x.units ** 0.5)


def fabs(x):
    x = _cq(x)
    return CQ(ca.fabs(x.expr) * x.units)


def sumsqr(x):
    x = _cq(x)
    return CQ(ca.sumsqr(x.expr) * x.units ** 2)


def sum1(x):
    x = _cq(x)
    return CQ(ca.sum1(x.expr) * x.units)


def mtimes(*args):
    r = args[0]
    for a in args[1:]:
        r = _cq(r) @ a
    return r


def _cat(cafun, args):
    args = [_cq(a) for a in args]
    u = args[0].units
    return CQ(cafun(*[a.m_as(u) for a in args]) * u)


def vertcat(*args):
    """Concatenate, converting everything to the unit of the first argument."""
    return _cat(ca.vertcat, args)


def horzcat(*args):
    return _cat(ca.horzcat, args)


def jacobian(y, x):
    """d y / d x, unit = unit(y) / unit(x).

    Only well defined when all entries of y (and of x) share one unit,
    which CQ enforces by construction."""
    y, x = _cq(y), _cq(x)
    return CQ(ca.jacobian(y.expr, x.expr) * (y.units / x.units))


def gradient(y, x):
    y, x = _cq(y), _cq(x)
    return CQ(ca.gradient(y.expr, x.expr) * (y.units / x.units))


def Function(name, ins, outs, *args, **kwargs):
    """ca.Function over CQ inputs/outputs; returns (Function, in_units, out_units).

    Inputs must be pure symbols, so they are passed in their own unit."""
    f = ca.Function(name, [_cq(i).expr for i in ins],
                    [_cq(o).expr for o in outs], *args, **kwargs)
    return f, [_cq(i).units for i in ins], [_cq(o).units for o in outs]


# Make pint hand mixed operations `Quantity <op> CQ` over to CQ (pint returns
# NotImplemented for registered "upcast" types, as it does for xarray), so
# e.g. `(A / ureg.h) @ x` yields a CQ instead of a Quantity wrapping a CQ.
try:
    from pint import compat as _pint_compat
    _pint_compat.upcast_type_map[_pint_compat.fully_qualified_name(CQ)] = CQ
except (ImportError, AttributeError):   # older/newer pint layout
    pass
