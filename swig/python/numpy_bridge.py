"""CasADi numpy NEP-13 (__array_ufunc__) and NEP-18 (__array_function__) bridge.

This module is spliced into the generated casadi.py via the
`%pythoncode "python/numpy_bridge.py"` directive in `swig/casadi.i`.
It runs in casadi module scope, so unqualified names like vcat/hcat/
if_else/floor/etc. resolve against casadis own exports, and `_casadi`
refers to the underlying C-extension submodule.

Extracted from casadi.i in commit ...; see commits d9e4ee1de1,
be1ad4c5a4, 2fd38a81a2, bef9994f81 for the four rounds of handler
authoring that built this bridge.

For new handler authoring, see the skills:
* numpy_nep13_18_kwargs_hygiene -- signature/kwargs discipline
* casadi_numpy_bridge_handler_authoring -- vocabulary traps
"""

# ------------------------------------------------------------------
# Block 1: NEP-13 ufunc dispatch (__array_ufunc__)
# ------------------------------------------------------------------
import math

try:
  from numpy import pi, inf, sum
except:
  pass

_LN2 = math.log(2.0)
_RAD_PER_DEG = math.pi / 180.0
_DEG_PER_RAD = 180.0 / math.pi

arcsin = lambda x: _casadi.asin(x)
arccos = lambda x: _casadi.acos(x)
arctan = lambda x: _casadi.atan(x)
arctan2 = lambda x,y: _casadi.atan2(x, y)
arctanh = lambda x: _casadi.atanh(x)
arcsinh = lambda x: _casadi.asinh(x)
arccosh = lambda x: _casadi.acosh(x)


def _pairwise_reduce(op, items):
    """Reduce `items` with the binary `op` using a balanced tree.

    Casadi has no native `prod` / per-axis `mmax` reduction, so handlers
    that need them fold the operation in Python.  A naive left-deep
    chain (`r = items[0]; for x in items[1:]: r = op(r, x)`) builds a
    symbolic graph of depth O(N), which produces deep gradients and
    slows down both AD and Function construction.  A balanced tree has
    depth O(log N) with the same total node count -- strictly better.

    Returns a single value of the same type as the inputs.  Caller
    must ensure `items` is non-empty.
    """
    items = list(items)
    if not items:
        raise ValueError("_pairwise_reduce: empty input")
    while len(items) > 1:
        n = len(items)
        out = [op(items[i], items[i + 1]) for i in range(0, n - 1, 2)]
        if n & 1:
            out.append(items[-1])
        items = out
    return items[0]


def _np_modf_ufunc(x):
    # numpy.modf returns (fractional, integer) parts, with integer-part
    # truncation toward zero (sign(x) * floor(|x|)).  if_else gives a
    # branchless symbolic form.
    int_part = _casadi.if_else(x >= 0, _casadi.floor(x), _casadi.ceil(x))
    return (x - int_part, int_part)


def _np_logaddexp2_ufunc(x, y):
    # log2(2**x + 2**y) with the standard max-shift for stability.
    m = _casadi.fmax(x, y)
    delta = _casadi.fmin(x, y) - m
    return m + _casadi.log1p(_casadi.exp(delta * _LN2)) / _LN2


def _np_logaddexp_ufunc(x, y):
    # log(exp(x) + exp(y)) -- delegate to casadi.logsumexp(vertcat(x,y))
    # for scalar inputs so MX gets the single OP_LOGSUMEXP node (with
    # the proven stable max-shift and closed-form softmax derivative).
    # For element-wise vector / matrix inputs casadi.logsumexp would
    # reduce the whole stack to a scalar, so fall back to the
    # element-wise max-shift formulation.
    def _is_scalar(v):
        return not hasattr(v, "shape") or v.shape == (1, 1)
    if _is_scalar(x) and _is_scalar(y):
        return _casadi.logsumexp(vertcat(x, y))
    m = _casadi.fmax(x, y)
    return m + _casadi.log1p(_casadi.exp(_casadi.fmin(x, y) - m))


def _np_divmod_ufunc(x, y):
    # numpy.divmod returns (floor_div, mod) -- pythonic floor-division,
    # not C-style truncation; mod has the sign of y.
    q = _casadi.floor(x / y)
    return (q, x - y * q)


def _np_heaviside_ufunc(x, h):
    # numpy.heaviside(x, h) = 0 if x<0, h if x==0, 1 if x>0.
    return _casadi.if_else(x > 0, 1.0, _casadi.if_else(x < 0, 0.0, h))


# numpy ufunc.__name__  ->  callable mapping casadi/numeric inputs to
# a casadi result.  Single source of truth; adding a new ufunc means
# adding one line.  Used by DM/SX/MX.__array_ufunc__ via
# _numpy_ufunc_dispatch below.
_NUMPY_UFUNC_DISPATCH = {
    # --- unary ---
    "negative":      lambda x: -x,
    "positive":      lambda x:  x,
    "absolute":      _casadi.fabs,
    "fabs":          _casadi.fabs,
    "sqrt":          _casadi.sqrt,
    "square":        lambda x: x*x,
    "reciprocal":    lambda x: 1/x,
    "sign":          _casadi.sign,
    "exp":           _casadi.exp,
    "expm1":         _casadi.expm1,
    "log":           _casadi.log,
    "log1p":         _casadi.log1p,
    "log10":         _casadi.log10,
    "sin":           _casadi.sin,
    "cos":           _casadi.cos,
    "tan":           _casadi.tan,
    "arcsin":        _casadi.asin,
    "arccos":        _casadi.acos,
    "arctan":        _casadi.atan,
    "sinh":          _casadi.sinh,
    "cosh":          _casadi.cosh,
    "tanh":          _casadi.tanh,
    "arcsinh":       _casadi.asinh,
    "arccosh":       _casadi.acosh,
    "arctanh":       _casadi.atanh,
    "floor":         _casadi.floor,
    "ceil":          _casadi.ceil,
    "logical_not":   _casadi.logic_not,
    # --- binary ---
    "add":           _casadi.plus,
    "subtract":      _casadi.minus,
    "multiply":      _casadi.times,
    "divide":        _casadi.rdivide,
    "true_divide":   _casadi.rdivide,
    "power":         _casadi.power,
    "float_power":   _casadi.power,
    "fmin":          _casadi.fmin,
    "fmax":          _casadi.fmax,
    "minimum":       _casadi.fmin,
    "maximum":       _casadi.fmax,
    # numpy.fmod follows C: sign of dividend.  Maps to casadi.fmod.
    "fmod":          _casadi.fmod,
    # numpy.mod / numpy.remainder follow Python %: sign of divisor.
    # casadi.remainder is IEEE-754 (round-to-nearest); we need a-b*floor(a/b).
    "remainder":     lambda x, y: x - y * _casadi.floor(x / y),
    "arctan2":       _casadi.atan2,
    "hypot":         _casadi.hypot,
    "copysign":      _casadi.copysign,
    "less":          _casadi.lt,
    "less_equal":    _casadi.le,
    "greater":       lambda x, y: _casadi.lt(y, x),
    "greater_equal": lambda x, y: _casadi.le(y, x),
    "equal":         _casadi.eq,
    "not_equal":     _casadi.ne,
    "logical_and":   _casadi.logic_and,
    "logical_or":    _casadi.logic_or,
    "matmul":        _casadi.mtimes,
    # --- batch 2: synthesized from casadi primitives (no native equivalent) ---
    "cbrt":          lambda x: _casadi.sign(x) * _casadi.power(_casadi.fabs(x), 1.0/3.0),
    "exp2":          lambda x: _casadi.exp(x * _LN2),
    "log2":          lambda x: _casadi.log(x) / _LN2,
    "trunc":         lambda x: _casadi.sign(x) * _casadi.floor(_casadi.fabs(x)),
    # NB: rint() is round-half-AWAY-from-zero, NOT numpy's banker's rounding.
    # Tie-breaking differs at exact half-integers; all other values match.
    "rint":          lambda x: _casadi.sign(x) * _casadi.floor(_casadi.fabs(x) + 0.5),
    "deg2rad":       lambda x: x * _RAD_PER_DEG,
    "radians":       lambda x: x * _RAD_PER_DEG,
    "degrees":       lambda x: x * _DEG_PER_RAD,
    "rad2deg":       lambda x: x * _DEG_PER_RAD,
    # signbit() returns 0/1 DM (not bool) for DM inputs; close enough for
    # numpy-compat under np.allclose().  -0.0 vs +0.0 distinction is lost.
    "signbit":       lambda x: _casadi.lt(x, 0),
    "floor_divide":  lambda x, y: _casadi.floor(x / y),
    "modf":          _np_modf_ufunc,
    "logaddexp2":    _np_logaddexp2_ufunc,
    # --- batch 3: round 3 ufuncs ---
    "logaddexp":     _np_logaddexp_ufunc,
    "divmod":        _np_divmod_ufunc,
    "heaviside":     _np_heaviside_ufunc,
    # ldexp(x, n) = x * 2**n; n is integer in numpy, but symbolic n is OK in casadi.
    "ldexp":         lambda x, n: x * _casadi.power(2.0, n),
    # --- batch 4: round 4 ufuncs ---
    # logical_xor: casadi has no native; coerce both operands to 0/1 then compare.
    "logical_xor":   lambda a, b: _casadi.ne(_casadi.ne(a, 0), _casadi.ne(b, 0)),
    # NB: isinf / isnan / isfinite are intentionally NOT registered.
    # For DM the missing-entry path falls back to numpy on densified
    # inputs (correct).  For SX/MX a "return symbolic 0" handler would
    # be unsound -- the expression might evaluate to NaN/inf at runtime
    # -- so we let the fallback's TypeError propagate as a clear error.
}

def _np_reduce_sum(x, axis):
    if axis is None:
        return _casadi.sum(x)
    if axis == 0:
        return _casadi.sum1(x)
    if axis == 1:
        return _casadi.sum2(x)
    raise ValueError("axis {0} is out of bounds for casadi 2-D array".format(axis))

def _np_reduce_global(fn):
    def go(x, axis):
        if axis is not None:
            raise NotImplementedError(
                "casadi {0} only supports an axis-less reduction".format(fn.__name__))
        return fn(x)
    return go


def _np_reduce_logaddexp(x, axis):
    # np.logaddexp.reduce -> log-sum-exp.  casadi has a stable, AD-friendly
    # logsumexp(column-vector).  For axis=None on a 2-D input we flatten
    # to a column first; for axis=0/1 we synthesize log(sum_axis(exp(...)))
    # with a max-shift so the result stays AD-friendly and numerically
    # well-behaved.
    if axis is None:
        flat = _casadi.vec(_casadi.transpose(x)) \
                if x.shape[1] != 1 else _casadi.vec(x)
        return _casadi.logsumexp(flat)
    if axis in (0, -2):
        m = _casadi.mmax(x)  # global max (conservative, still stable)
        return m + _casadi.log(_casadi.sum1(_casadi.exp(x - m)))
    if axis in (1, -1):
        m = _casadi.mmax(x)
        return m + _casadi.log(_casadi.sum2(_casadi.exp(x - m)))
    raise NotImplementedError(
        "axis {0} not supported for logaddexp.reduce".format(axis))


# numpy ufunc.__name__  ->  (reduce_handler, accumulate_handler) where each
# takes (casadi_value, axis) and returns a casadi result (or raises to fall
# back to numpy / give a clear error).
_NUMPY_UFUNC_REDUCE = {
    "add":         (_np_reduce_sum,
                    lambda x, axis: _casadi.cumsum(x, 0 if axis is None else axis)),
    "maximum":     (_np_reduce_global(_casadi.mmax), None),
    "minimum":     (_np_reduce_global(_casadi.mmin), None),
    "fmax":        (_np_reduce_global(_casadi.mmax), None),
    "fmin":        (_np_reduce_global(_casadi.mmin), None),
    "logical_and": (_np_reduce_global(_casadi.logic_all), None),
    "logical_or":  (_np_reduce_global(_casadi.logic_any), None),
    # numpy: np.logaddexp.reduce(x) == log(sum(exp(x))); casadi has a
    # native, stable logsumexp(column-vector) and we synthesize the
    # axis-aware variants on top.  Wires through to symbolic SX/MX as
    # well -- the fallback densifier can only handle DM.
    "logaddexp":   (_np_reduce_logaddexp, None),
}


def _numpy_ufunc_dispatch(self, ufunc, method, inputs, kwargs):
    # NEP-13 __array_ufunc__ entry point, shared by DM/SX/MX.
    name = ufunc.__name__
    if name.endswith(" (vectorized)"):       # numpy.frompyfunc decoration
        name = name[:-len(" (vectorized)")]

    # Gateway (issue #2959).  mode 1: the casadi-aware numpy support routes
    # the whole call through the numpy-semantics wrapper, which returns a
    # NumpyArray following numpy's shape/axis contract -- for SX/MX/DM alike.
    # mode 0/-1: legacy casadi 3.7.2 behaviour (0 warns, -1 is silent).
    _mode = GlobalOptions.getNumpyMode()
    if _mode == 1:
        return ArrayInterface._wrap(self, 2).__array_ufunc__(ufunc, method, *inputs, **kwargs)

    # Default (legacy casadi 3.7.2) mode.
    # Operator ufuncs called as an operator (`x + M`, `-M`, `np_arr < M`)
    # keep returning casadi types silently: ordinary operator interop, not
    # the explicit `numpy.foo(M)` surface the deprecation targets.
    if method == "__call__" and name in _OPERATOR_UFUNCS:
        op = _NUMPY_UFUNC_DISPATCH.get(name)
        if op is not None:
            try:
                return op(*inputs)
            except (TypeError, ValueError, NotImplementedError):
                pass
        return _numpy_ufunc_fallback(self, ufunc, method, inputs, kwargs)

    # Any other explicit numpy ufunc call (np.sin(M), np.add.reduce(M), ...).
    # 3.7.2 behaviour: numeric inputs densify to a numpy result; symbolic
    # inputs go through the casadi op and return a casadi value (np.sin(MX)
    # -> MX).  Mode 0 also emits a FutureWarning; mode -1 stays silent.
    if _mode == 0:
        _warn_preserve_type()
    if not _has_symbolic(inputs):
        for k in ("keepdims", "casting", "order", "dtype", "subok", "signature"):
            kwargs.pop(k, None)
        return _numpy_ufunc_fallback(self, ufunc, method, inputs, kwargs)
    # Symbolic: apply the casadi-native op / reduce (3.7.2 returned casadi).
    if method == "__call__":
        op = _NUMPY_UFUNC_DISPATCH.get(name)
        if op is not None:
            try:
                return op(*inputs)
            except (TypeError, ValueError, NotImplementedError):
                pass
    elif method in ("reduce", "accumulate"):
        entry = _NUMPY_UFUNC_REDUCE.get(name)
        idx = 0 if method == "reduce" else 1
        if entry is not None and entry[idx] is not None:
            try:
                return entry[idx](inputs[0], kwargs.get("axis", 0))
            except (TypeError, ValueError, NotImplementedError):
                pass
    return _numpy_ufunc_fallback(self, ufunc, method, inputs, kwargs)


def _numpy_ufunc_fallback(self, ufunc, method, inputs, kwargs):
    # DM: round-trip via .full() to numpy and re-dispatch.
    # Symbolic types (SX/MX): raise; there is no meaningful conversion.
    if not hasattr(self, "full"):
        raise TypeError(
            "casadi {0} does not support numpy.{1}.{2}.  "
            "Use the equivalent casadi function instead.".format(
                type(self).__name__, ufunc.__name__, method))
    new_inputs = tuple(x.full() if hasattr(x, "full") else x for x in inputs)
    return getattr(ufunc, method)(*new_inputs, **kwargs)

# ------------------------------------------------------------------
# Block 2: NEP-18 free-function dispatch (__array_function__)
# ------------------------------------------------------------------

# NEP 18 (__array_function__) bridge.
#
# Lets casadi values participate in calls like np.concatenate, np.where,
# np.reshape, np.linalg.solve, ...  numpy invokes __array_function__ on
# the first argument that defines it; we look the function up in
# _NUMPY_FUNCTION_DISPATCH below and either delegate to the casadi
# equivalent or return NotImplemented (so numpy falls back to its
# default).  The table is built lazily on first dispatch because the
# numpy module isn't guaranteed to be importable at SWIG-init time.

_NUMPY_FUNCTION_DISPATCH = None

def _is_casadi_value(x):
    return isinstance(x, (DM, SX, MX))

def _has_symbolic(values):
    # True if any operand is a symbolic casadi value (SX/MX).
    return any(isinstance(x, (SX, MX)) for x in values)

# Numpy ufuncs that correspond to a Python operator (`x + M`, `-M`,
# `np_arr < M`, ...).  When invoked as `__call__` in DEFAULT (legacy)
# mode these keep returning casadi types silently -- that is ordinary
# casadi/numpy operator interop, not the explicit `numpy.foo(M)` surface
# the issue #2959 deprecation targets.  (Reductions like `np.add.reduce`
# use method="reduce", so they fall through to the deprecation path.)
# Ufuncs reachable through a Python binary operator / comparison / builtin
# (`arr + M`, `arr < M`, `divmod(arr, M)`): there `np.foo(arr, M)` and the
# operator are the SAME __array_ufunc__ call, so these must stay silent or
# ordinary interop would warn.  The UNARY operator ufuncs are deliberately
# NOT here: `-M` / `abs(M)` go via __neg__ / __abs__ and never reach
# __array_ufunc__, so an explicit np.negative(M) / np.abs(M) is
# distinguishable and warns like np.sin.
_OPERATOR_UFUNCS = frozenset([
    "add", "subtract", "multiply", "matmul", "divide", "true_divide",
    "floor_divide", "remainder", "mod", "divmod", "power", "float_power",
    "less", "less_equal", "greater", "greater_equal", "equal", "not_equal",
])

# A single CONSTANT notice text, always emitted from the SAME line below.
# That makes the stdlib's default warning action dedupe it to once-per-
# process on its own -- no custom flag and, crucially, NO global warnings-
# filter mutation, so other packages' warnings are untouched.  The standard
# overrides still apply: `python -W always` / PYTHONWARNINGS=always shows
# every occurrence, `-W error` makes it fatal.  (We deliberately do NOT use a
# user-call-site stacklevel: a varying location would defeat the dedupe and
# flood multi-site scripts.)
_NUMPY_LEGACY_NOTICE = (
    "\n"
    "casadi: a numpy function was called on a casadi value (issue #2959).\n"
    "This used legacy casadi 3.7.2 behaviour, where the result type follows\n"
    "the input:\n"
    "  - a numeric input (DM)       -> a plain numpy array\n"
    "  - a symbolic input (SX / MX) -> a casadi value, or an error if the\n"
    "                                  function has no casadi equivalent\n"
    "\n"
    "For a shape-correct, type-preserving result, opt in to the\n"
    "casadi-aware numpy support:\n"
    "\n"
    "    try:\n"
    "        ca.GlobalOptions.setNumpyMode(1)\n"
    "    except AttributeError:\n"
    "        pass\n"
    "\n"
    "Or use setNumpyMode(-1) to keep legacy behaviour silently.\n"
    "Shown once; run python with -W always to see every occurrence.\n")


def _warn_preserve_type():
    # issue #2959: an explicit numpy.foo(M) on a casadi value uses legacy
    # 3.7.2 behaviour by default (numeric densifies, symbolic returns casadi).
    # FutureWarning (not Deprecation): it changes the RESULT, so per numpy
    # convention it must be visible by default; DeprecationWarning would be
    # suppressed outside __main__.  Default stacklevel -> this constant line,
    # which is what lets the stdlib dedupe to once-per-process.
    import warnings as _w
    _w.warn(_NUMPY_LEGACY_NOTICE, FutureWarning)


def _np_concatenate(arrays, axis=0):
    arrays = list(arrays)
    if axis == 0 or axis is None:
        return vcat(arrays)
    if axis == 1:
        return hcat(arrays)
    return NotImplemented

def _np_vstack(tup):
    return vcat(tup)

def _np_hstack(tup):
    return hcat(tup)

def _np_stack(arrays, axis=0):
    if axis == 0:
        return vcat([a.T if hasattr(a, "T") else a for a in arrays])
    if axis == 1:
        return hcat(arrays)
    return NotImplemented

def _np_column_stack(tup):
    return hcat([a if hasattr(a, "shape") and len(a.shape) > 1 and a.shape[1] > 0 else
                     reshape(a, (-1, 1)) for a in tup])

def _np_where(condition, x=None, y=None):
    if x is None and y is None:
        return NotImplemented   # nz-index variant: leave to numpy
    # short_circuit=False -> elementwise (cond*x + (1-cond)*y) semantics,
    # which is what numpy.where does.  short_circuit=True builds a Switch
    # node that requires a scalar condition.
    return if_else(condition, x, y, False)

def _np_reshape(a, newshape, order='C'):
    # casadi.reshape is column-major (Fortran).  numpy default is C-order
    # (row-major).  Transpose-reshape-transpose to get numpy semantics.
    if isinstance(newshape, int):
        newshape = (newshape, 1)
    if len(newshape) == 1:
        newshape = (newshape[0], 1)
    m, n = int(newshape[0]), int(newshape[1])
    if order == 'F':
        return reshape(a, m, n)
    if order == 'C':
        return reshape(a.T, n, m).T
    return NotImplemented

def _np_transpose(a, axes=None):
    if axes is not None and tuple(axes) not in ((0, 1), (1, 0)):
        return NotImplemented
    return a.T

def _np_dot(a, b, out=None):
    if out is not None:
        return NotImplemented
    if hasattr(a, "is_vector") and a.is_vector() and hasattr(b, "is_vector") and b.is_vector():
        return dot(a, b)
    return mtimes(a, b)

def _np_inner(a, b):
    return dot(a, b)


def _np_einsum(subscripts, *operands, out=None, dtype=None, order='K',
               casting='safe', optimize=False):
    # Bridge to casadi.einstein, which takes flat tensors plus per-
    # label dim lists and integer label codes.  Supports the common
    # 2-operand cases (matmul, mat-vec, vec-mat, inner product,
    # outer-product-via-explicit-axes, transpose-of-product).
    # Symbolic SX/MX einsum was previously a hard "no implementation
    # found" -- this handler lights up that path.
    if out is not None or dtype is not None:
        return NotImplemented
    if not isinstance(subscripts, str):
        return NotImplemented   # interleaved-list form
    if '->' not in subscripts:
        return NotImplemented   # implicit output requires repeated-
                                # label inference -- caller can use
                                # explicit "->" or fall back to numpy.
    in_str, out_str = subscripts.split('->', 1)
    in_lbls = in_str.split(',')
    if len(in_lbls) != len(operands):
        return NotImplemented
    if len(operands) != 2:
        return NotImplemented   # casadi.einstein is binary/ternary;
                                # unary einsum (transpose, diagonal,
                                # axis-reduce) doesn't map to it.
    if not all(hasattr(op, 'shape') for op in operands):
        return NotImplemented

    def _per_label_dims(op, lbls):
        # 2 labels -> use shape as-is (casadi values are always 2-D).
        # 1 label  -> require a vector shape (one dim == 1); collapse.
        # 0 labels -> scalar input, not bridgeable through einstein.
        n = len(lbls)
        if n == 2:
            return list(op.shape)
        if n == 1:
            r, c = op.shape
            if r == 1:
                return [c]
            if c == 1:
                return [r]
            return None
        return None

    da = _per_label_dims(operands[0], in_lbls[0])
    db = _per_label_dims(operands[1], in_lbls[1])
    if da is None or db is None:
        return NotImplemented

    # Stable letter -> negative-integer label mapping.
    seen = []
    for ch in in_str.replace(",", "") + out_str:
        if ch not in seen:
            seen.append(ch)
    letter_to_label = {ch: -(i + 1) for i, ch in enumerate(seen)}

    a_lbls = [letter_to_label[ch] for ch in in_lbls[0]]
    b_lbls = [letter_to_label[ch] for ch in in_lbls[1]]
    c_lbls = [letter_to_label[ch] for ch in out_str]

    # Per-letter dim: consistent across inputs.
    letter_dim = {}
    for lbls, dims in ((in_lbls[0], da), (in_lbls[1], db)):
        for letter, d in zip(lbls, dims):
            if letter in letter_dim and letter_dim[letter] != d:
                return NotImplemented
            letter_dim[letter] = d
    if not all(ch in letter_dim for ch in out_str):
        return NotImplemented
    dc = [letter_dim[ch] for ch in out_str]

    Av = vec(operands[0])
    Bv = vec(operands[1])
    flat = _casadi.einstein(Av, Bv, da, db, dc, a_lbls, b_lbls, c_lbls)
    if len(dc) == 0:
        return flat                    # scalar
    if len(dc) == 1:
        return reshape(flat, dc[0], 1) # column-vector shape
    if len(dc) == 2:
        return reshape(flat, dc[0], dc[1])
    return NotImplemented   # >2 output labels not representable in 2-D


def _np_outer(a, b, out=None):
    if out is not None:
        return NotImplemented
    return mtimes(reshape(a, (-1, 1)), reshape(b, (1, -1)))

def _np_cumsum(a, axis=None):
    return cumsum(a, 0 if axis is None else int(axis))

def _np_repeat(a, repeats, axis=None):
    # numpy.repeat repeats each element along an axis (or flattens, in C
    # order, when axis is None).  Equivalent to a Kronecker product with
    # an all-ones vector of the appropriate shape.
    if not isinstance(repeats, int):
        return NotImplemented   # numpy supports per-element repeat counts
    n = int(repeats)
    # Result orientation: numpy returns 1-D for axis=None.  We follow the
    # casadi convention "1-D becomes column" for the axis=None case.
    if axis is None:
        return kron(vec(a.T), DM.ones(n, 1))     # row-major flatten then repeat
    if axis == 0:
        return kron(a, DM.ones(n, 1))
    if axis == 1:
        return kron(a, DM.ones(1, n))
    return NotImplemented

def _np_tile(a, reps):
    # numpy.tile semantics for a 2-D input: an integer or a length-1
    # tuple tiles along the *last* axis.  A length-2 tuple tiles
    # (rows, cols).  casadi DMs are always 2-D, so the 2-D case applies.
    if isinstance(reps, int):
        return repmat(a, 1, int(reps))
    if len(reps) == 1:
        return repmat(a, 1, int(reps[0]))
    if len(reps) == 2:
        return repmat(a, int(reps[0]), int(reps[1]))
    return NotImplemented

def _np_zeros_like(a):
    return type(a).zeros(a.sparsity()) if hasattr(a, "sparsity") else NotImplemented

def _np_ones_like(a):
    return type(a).ones(a.sparsity()) if hasattr(a, "sparsity") else NotImplemented

def _np_full_like(a, fill_value):
    return type(a)(a.sparsity(), fill_value) if hasattr(a, "sparsity") else NotImplemented

def _np_sum(a, axis=None, keepdims=False):
    # keepdims is accepted but ignored: casadi values are always 2-D, so
    # the result already matches numpy's keepdims=True shape.
    if axis is None:
        return sum1(sum2(a))     # total sum, kept symbolic
    if axis == 0:
        return sum1(a)
    if axis == 1:
        return sum2(a)
    raise ValueError("axis {0} is out of bounds for casadi 2-D array".format(axis))

def _np_max(a, axis=None, keepdims=False):
    if axis is None:
        return mmax(a)
    return NotImplemented

def _np_min(a, axis=None, keepdims=False):
    if axis is None:
        return mmin(a)
    return NotImplemented

def _np_all(a, axis=None, keepdims=False):
    if axis is None:
        return logic_all(a)
    return NotImplemented

def _np_any(a, axis=None, keepdims=False):
    if axis is None:
        return logic_any(a)
    return NotImplemented

def _np_linalg_norm(a, ord=None, axis=None, keepdims=False):
    # Matches numpy.linalg.norm semantics on 2-D inputs (which casadi
    # values always are).  In numpy this means:
    #   ord=None / 'fro'  ->  Frobenius
    #   ord=1             ->  max column sum
    #   ord=-1            ->  min column sum
    #   ord=inf           ->  max row sum
    #   ord=-inf          ->  min row sum
    #   ord=2             ->  spectral norm (largest singular value)
    # casadi's norm_1 / norm_inf are element-wise (sum-of-abs /
    # max-of-abs over all entries), NOT the induced matrix norms, so we
    # build the matrix norms from sum1/sum2/fabs/mmax/mmin.  For ord=2
    # there is no general casadi spectral norm, but for shape (n,1) or
    # (1,n) the spectral norm equals the vector 2-norm and casadi's
    # norm_2 handles those.
    if axis is not None or keepdims:
        return NotImplemented
    inf = float("inf")
    is_vec_shape = hasattr(a, "is_vector") and a.is_vector()
    if ord is None or ord == "fro":
        return norm_fro(a)
    if ord == 2:
        return norm_2(a) if is_vec_shape else NotImplemented
    if ord == 1:
        return mmax(sum1(fabs(a)))
    if ord == -1:
        return mmin(sum1(fabs(a)))
    if ord in (inf, "inf"):
        return mmax(sum2(fabs(a)))
    if ord in (-inf, "-inf"):
        return mmin(sum2(fabs(a)))
    return NotImplemented


def _np_clip(a, a_min=None, a_max=None):
    if a_min is not None:
        a = fmax(a, a_min)
    if a_max is not None:
        a = fmin(a, a_max)
    return a


def _np_diff(a, n=1, axis=-1):
    # casadi.diff(x, n, axis) supports repeated differencing and either
    # axis natively.  numpy's axis=-1 default means "last axis", which on
    # a casadi (n,1) column-vector would be the size-1 axis -> empty
    # (n,0).  Since casadi values are always 2-D, treat (n,1) and (1,n)
    # vector-shaped inputs as 1-D-like under the default axis: diff
    # along the populated axis instead of the size-1 one.  Explicit
    # axis= still honored literally.
    if not isinstance(n, int) or n < 0:
        return NotImplemented
    if axis == -1:
        if hasattr(a, "shape"):
            nr, nc = a.shape
            if nc == 1 and nr >= 1:
                axis = 0
            else:
                axis = 1
        else:
            axis = 1
    if axis not in (0, 1):
        return NotImplemented
    if n == 0:
        return a
    return diff(a, int(n), int(axis))


def _np_trapz(y, x=None, dx=1.0, axis=-1):
    """numpy.trapz handler: trapezoidal integration along an axis.

    On a casadi (n,1) or (1,n) value under the default axis=-1, treats
    the input as 1-D and integrates along the populated axis (matching
    the friendly convention used by `_np_diff` / `_np_gradient`).
    Returns a scalar for vector inputs; for matrices returns a vector
    of per-axis integrals.

    Note: numpy and aerosandbox.numpy disagree on what `trapz` means.
    Numpy's trapz returns the SUM of trapezoidal contributions (total
    integral); aerosandbox's trapz returns the per-subinterval values
    without summing.  We match numpy.
    """
    if not hasattr(y, "shape"):
        return NotImplemented
    nr, nc = y.shape
    is_col = (nc == 1 and nr >= 2)
    is_row = (nr == 1 and nc >= 2)
    if axis == -1:
        axis = 0 if is_col else 1
    if axis not in (0, 1):
        return NotImplemented
    n_along = nr if axis == 0 else nc
    if n_along < 2:
        return NotImplemented

    def take(a, lo, hi):
        if axis == 0:
            return a[lo:hi, :]
        return a[:, lo:hi]

    pair_sum = (take(y, 0, -1) + take(y, 1, None)) * 0.5
    if x is None:
        contributions = pair_sum * dx
    else:
        x_dm = x if hasattr(x, "shape") else DM(x)
        if x_dm.numel() != n_along:
            return NotImplemented
        x_v = vec(x_dm) if axis == 0 else vec(x_dm).T
        contributions = pair_sum * diff(x_v, 1, axis)
    return sum1(contributions) if axis == 0 else sum2(contributions)


def _np_gradient(f, *varargs, axis=None, edge_order=1):
    """numpy.gradient handler.

    Second-order accurate central differences in the interior;
    first-order one-sided differences at the boundaries (edge_order=1)
    or second-order one-sided at the boundaries (edge_order=2).  Treats
    a casadi (n,1) or (1,n) value as 1-D when axis is unspecified --
    same friendly convention as `_np_diff`, since casadi has no native
    1-D shape.  Spacing varargs accept either a scalar (uniform) or a
    1-D coordinate array of the right length.
    """
    if not hasattr(f, "shape"):
        return NotImplemented
    if edge_order not in (1, 2):
        return NotImplemented
    nr, nc = f.shape
    is_col = (nc == 1 and nr >= 2)
    is_row = (nr == 1 and nc >= 2)
    if axis is None:
        if is_col:
            axis = 0
        elif is_row:
            axis = 1
        else:
            # Matrix without explicit axis: numpy returns one gradient
            # per axis as a list.
            return [_np_gradient(f, *varargs, axis=ax, edge_order=edge_order)
                    for ax in (0, 1)]
    if isinstance(axis, (tuple, list)):
        return [_np_gradient(f, *varargs, axis=int(ax), edge_order=edge_order)
                for ax in axis]
    if axis not in (0, 1):
        return NotImplemented
    n_along = nr if axis == 0 else nc
    if n_along < 2:
        return NotImplemented
    if edge_order == 2 and n_along < 3:
        return NotImplemented

    # Resolve spacing along this axis.
    if len(varargs) == 0:
        h_uniform = 1.0
        coords = None
    elif len(varargs) == 1:
        s = varargs[0]
        s_dm = s if hasattr(s, "shape") else DM(s)
        if s_dm.numel() == 1:
            h_uniform = s_dm
            coords = None
        elif s_dm.numel() == n_along:
            h_uniform = None
            coords = vec(s_dm) if axis == 0 else vec(s_dm).T
        else:
            return NotImplemented
    else:
        # Per-axis varargs are a numpy multi-axis convenience; the
        # explicit-axis path above already picked one axis.
        return NotImplemented

    def take(a, lo, hi):
        if axis == 0:
            return a[lo:hi, :]
        return a[:, lo:hi]

    # Build interior spacings hm[i] (between i-1 and i) and hp[i]
    # (between i and i+1) for i in [1, n_along-2].
    if h_uniform is not None:
        hm = h_uniform
        hp = h_uniform
        hm_first = h_uniform
        hp_last = h_uniform
    else:
        dx_full = diff(coords, 1, axis)   # length n_along-1 along axis
        hm = take(dx_full, 0, -1)         # length n_along-2
        hp = take(dx_full, 1, None)       # length n_along-2
        hm_first = take(dx_full, 0, 1)
        hp_last = take(dx_full, -1, None)

    # Interior: 2nd-order central
    dfp = take(f, 2, None) - take(f, 1, -1)
    dfm = take(f, 1, -1)   - take(f, 0, -2)
    grad_int = (hm * hm * dfp + hp * hp * dfm) / (hm * hp * (hm + hp))

    # Edge gradients.
    if edge_order == 1:
        grad_first = (take(f, 1, 2) - take(f, 0, 1)) / hm_first
        grad_last  = (take(f, -1, None) - take(f, -2, -1)) / hp_last
    else:  # edge_order == 2
        f0 = take(f, 0, 1); f1 = take(f, 1, 2); f2 = take(f, 2, 3)
        fn = take(f, -1, None); fn1 = take(f, -2, -1); fn2 = take(f, -3, -2)
        if h_uniform is not None:
            h = h_uniform
            grad_first = (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h)
            grad_last  = ( 3.0 * fn - 4.0 * fn1 + fn2) / (2.0 * h)
        else:
            # Non-uniform 2nd-order at the boundaries; same closed form
            # aerosandbox uses, matches numpy's documented formula.
            dfm_f = f1 - f0; dfp_f = f2 - f1
            dfm_l = fn1 - fn2; dfp_l = fn - fn1
            hp_first = take(dx_full, 1, 2)       # between f[1] and f[2]
            hm_last  = take(dx_full, -2, -1)     # between f[-3] and f[-2]
            grad_first = (2.0 * dfm_f * hm_first * hp_first
                          + dfm_f * hp_first * hp_first
                          - dfp_f * hm_first * hm_first) \
                         / (hm_first * hp_first * (hm_first + hp_first))
            grad_last  = (-dfm_l * hp_last * hp_last
                          + dfp_l * hm_last * hm_last
                          + 2.0 * dfp_l * hm_last * hp_last) \
                         / (hm_last * hp_last * (hm_last + hp_last))

    if axis == 0:
        return vertcat(grad_first, grad_int, grad_last)
    return horzcat(grad_first, grad_int, grad_last)


def _np_diag(v, k=0):
    # numpy.diag has two modes (selected by input rank):
    #   1-D input  -> diagonal matrix with v on the k-th diagonal
    #                 (result shape (n+|k|, n+|k|))
    #   2-D input  -> extract the k-th diagonal as a 1-D array
    # CasADi values are always 2-D; we treat (n,1) and (1,n) as the
    # vector case (matching the convention used by other vec-like
    # handlers), and (m,n) with m,n>=2 as the matrix case.
    if not hasattr(v, "shape"):
        return NotImplemented
    is_vec = (v.shape[0] == 1 or v.shape[1] == 1)
    if is_vec:
        if k == 0:
            return diag(v)
        # Place v on the k-th diagonal of an (n+|k|, n+|k|) zero matrix.
        flat = vec(v) if v.shape[1] == 1 else vec(v.T)
        n = flat.numel()
        size = n + abs(k)
        base = diag(flat)                       # (n,n)
        if k > 0:
            return vertcat(horzcat(DM(n, k), base),
                           DM(k, size))
        # k < 0: shift down by |k| rows, pad right by |k| cols.
        return vertcat(DM(abs(k), size),
                       horzcat(base, DM(n, abs(k))))
    # Matrix case: extract k-th diagonal.
    m, n = v.shape
    if k == 0:
        # casadi.diag(matrix) extracts the main diagonal in one op for
        # square matrices; fall back to per-element vcat only when
        # non-square.
        if m == n:
            return diag(v)
        return vcat([v[i, i] for i in range(min(m, n))])
    if k > 0:
        # super-diagonals: v[i, i+k] for i in min(m, n-k).
        if k >= n:
            return DM.zeros(0, 1)
        return vcat([v[i, i + k] for i in range(min(m, n - k))])
    # k < 0: sub-diagonals: v[i-k, i] for i in min(n, m+k).
    if -k >= m:
        return DM.zeros(0, 1)
    return vcat([v[i - k, i] for i in range(min(n, m + k))])


def _np_split_offsets(n, indices_or_sections):
    # Translate numpy's `indices_or_sections` (int N or 1-D index iterable)
    # into a casadi-style offsets list spanning [0, n].  None signals
    # "unsupported" so the caller can propagate NotImplemented.
    import numpy as np
    if isinstance(indices_or_sections, (int, np.integer)):
        N = int(indices_or_sections)
        if N <= 0 or n % N != 0:
            return None
        return list(range(0, n + 1, n // N))
    try:
        idx = [int(i) for i in indices_or_sections]
    except (TypeError, ValueError):
        return None
    # numpy clamps split indices to [0, n] and inserts empty splits for
    # out-of-range; casadi tolerates that via repeated offsets.
    clipped = [min(max(i, 0), n) for i in idx]
    # casadi requires non-decreasing offsets; numpy lets indices be
    # unsorted but only the sorted case has well-defined semantics.
    if any(clipped[i] > clipped[i + 1] for i in range(len(clipped) - 1)):
        return None
    return [0] + clipped + [n]


def _np_hsplit(ary, indices_or_sections):
    offsets = _np_split_offsets(ary.shape[1], indices_or_sections)
    if offsets is None:
        return NotImplemented
    return horzsplit(ary, offsets)


def _np_vsplit(ary, indices_or_sections):
    offsets = _np_split_offsets(ary.shape[0], indices_or_sections)
    if offsets is None:
        return NotImplemented
    return vertsplit(ary, offsets)


def _np_trace(a, offset=0, axis1=0, axis2=1, dtype=None, out=None):
    # casadi.trace is the main-diagonal sum of a 2-D matrix.  numpy.trace
    # supports an offset and arbitrary axis pair; the latter is meaningless
    # for casadi's 2-D-only values, but the former is a real feature gap.
    if offset != 0 or out is not None or dtype is not None:
        return NotImplemented
    if (axis1, axis2) not in ((0, 1), (1, 0)):
        return NotImplemented
    if not hasattr(a, "shape"):
        return NotImplemented
    return trace(a)


def _np_linspace(start, stop, num=50, endpoint=True, retstep=False,
                 dtype=None, axis=0):
    # casadi.linspace generates `num` points from start to stop INCLUSIVE.
    # We extend to numpy's endpoint=False by generating num+1 points and
    # dropping the last.  retstep is honoured symbolically -- routing to
    # the numpy fallback with casadi-native endpoints produces a 3-D
    # shape (the (1,1)*(1,1) broadcast trap) so we must not delegate.
    if dtype is not None or axis != 0:
        return NotImplemented
    n = int(num)
    if n <= 0:
        return NotImplemented
    if endpoint:
        samples = linspace(start, stop, n)
        step = (stop - start) / (n - 1) if n > 1 else (stop - start)
    else:
        full = linspace(start, stop, n + 1)
        samples = full[:n]
        step = (stop - start) / n
    if retstep:
        return (samples, step)
    return samples


def _np_cross(a, b, axisa=-1, axisb=-1, axisc=-1, axis=None):
    # casadi.cross supports 3-component vectors only and takes a 1-indexed
    # `dim` argument (1 = vectors along rows, 2 = along columns, -1 = auto).
    # numpy.cross uses 0-indexed axes; `axis` overrides axisa/axisb/axisc.
    if axis is not None:
        axisa = axisb = axisc = axis
    if axisa != axisb or axisa != axisc:
        # casadi has no axis-permutation step; refuse mixed axes.
        return NotImplemented
    ax = axisa
    if ax in (-1, 1):
        dim = 2
    elif ax == 0:
        dim = 1
    else:
        return NotImplemented
    # ax=-1 is numpy default: prefer auto-detect so 1-D/column-vector
    # inputs (which casadi represents as (3,1)) keep working.
    if axis is None and axisa == -1:
        dim = -1
    return cross(a, b, dim)


def _np_roll(a, shift, axis=None):
    # 1-D roll equivalent: vec the input, slice and concat.
    if axis is None:
        flat = vec(a)
        n = flat.shape[0]
        s = int(shift) % n if n else 0
        if s == 0:
            return reshape(flat, a.shape[0], a.shape[1])
        rolled = vertcat(flat[n - s:], flat[:n - s])
        return reshape(rolled, a.shape[0], a.shape[1])
    if axis == 0:
        m = a.shape[0]
        if m == 0:
            return a
        s = int(shift) % m
        if s == 0:
            return a
        return vertcat(a[m - s:, :], a[:m - s, :])
    if axis == 1:
        n = a.shape[1]
        if n == 0:
            return a
        s = int(shift) % n
        if s == 0:
            return a
        return horzcat(a[:, n - s:], a[:, :n - s])
    return NotImplemented


def _np_atleast_1d(*arys):
    # casadi values are always >= 2-D; pass through unchanged, in line
    # with the NEP-18 convention used by dask, jax, xarray, PyData/sparse,
    # and scipy.sparse arrays.  Consumers that need the dense numpy form
    # (e.g. matplotlib) reach DM via .to_numpy() instead -- silent
    # densification here would defeat the whole point of preserving
    # sparse representations through shape ops.
    return arys[0] if len(arys) == 1 else list(arys)


def _np_atleast_2d(*arys):
    return arys[0] if len(arys) == 1 else list(arys)

# --- batch 2: synthesized handlers ---

def _np_arange(start, stop=None, step=1, dtype=None, device=None, like=None):
    import numpy as np
    # dtype: ignored — CasADi DM is float64-only.
    # device: ignored — CasADi has no device-placement concept.
    # like: ignored — protocol-only, has no semantic effect once we've
    #   already dispatched into this handler.
    if device is not None or like is not None:
        return NotImplemented

    # Normalize: numpy.arange(stop) form vs numpy.arange(start, stop[, step]).
    if stop is None:
        start, stop = 0, start

    # The count of elements depends on the numeric value of stop-start;
    # if any argument is a symbolic SX/MX, the size is non-deterministic.
    def _to_scalar(v):
        if hasattr(v, "is_constant"):
            if not v.is_constant():
                return None
            return float(v)
        return float(v)

    s = _to_scalar(start)
    e = _to_scalar(stop)
    d = _to_scalar(step)
    if s is None or e is None or d is None:
        return NotImplemented

    vals = np.arange(s, e, d)
    return DM(vals.reshape(-1, 1))

def _np_asanyarray(a, dtype=None, order=None, like=None):
    import numpy as np
    # dtype:  if a non-float dtype is requested for a numeric DM, hand off.
    # order:  C/F/A/K ordering -- meaningless for casadi; accept defaulted.
    # like:   array-API protocol arg; only None makes sense.
    if like is not None:
        return NotImplemented
    if dtype is not None:
        try:
            np_dtype = np.dtype(dtype)
        except TypeError:
            return NotImplemented
        if not np.issubdtype(np_dtype, np.floating):
            return NotImplemented
    return a

def _np_asarray(a, dtype=None, order=None, copy=None, like=None):
    import numpy as np
    # dtype:  cast request; only float-like accepted symbolically.
    # order:  C/F/A/K -- meaningless for casadi 2-D values.
    # copy:   True/False/None -- casadi values are reference-semantics;
    #         we ignore (the user gets back the same value).
    # like:   array-API protocol arg; only None.
    if like is not None:
        return NotImplemented
    if dtype is not None:
        try:
            np_dtype = np.dtype(dtype)
        except TypeError:
            return NotImplemented
        if not np.issubdtype(np_dtype, np.floating):
            return NotImplemented
    return a

def _np_ascontiguousarray(a, dtype=None, like=None):
    import numpy as np
    # dtype: metadata-only for symbolic types; if a non-float-like
    #        dtype is requested, hand off to fallback (DM can densify).
    # like:  array-API protocol arg; only None makes sense here.
    if like is not None:
        return NotImplemented
    if dtype is not None:
        try:
            np_dtype = np.dtype(dtype)
        except TypeError:
            return NotImplemented
        if not np.issubdtype(np_dtype, np.floating):
            return NotImplemented
    return a

def _np_bmat(obj, ldict=None, gdict=None):
    # ldict/gdict: only used by the string-spec form; we don't support
    #              string specs symbolically.
    if isinstance(obj, str):
        return NotImplemented
    if ldict is not None or gdict is not None:
        return NotImplemented
    if not isinstance(obj, (list, tuple)):
        return NotImplemented
    # Two forms:
    #   1) list-of-lists: bmat([[A, B], [C, D]]) -> vcat row-stacks
    #   2) flat list: bmat([A, B]) -> hcat (numpy treats as 1-D)
    if all(isinstance(row, (list, tuple)) for row in obj):
        rows = [hcat(list(row)) for row in obj]
        return vcat(rows)
    return hcat(list(obj))

def _np_copy(a, order='K', subok=False):
    # order: no-op. CasADi has a single storage convention (CCS), so the
    #   K/A/C/F memory-layout hints don't apply.
    # subok: no-op. CasADi DM/SX/MX have no ndarray-subclass hierarchy
    #   for numpy to preserve.
    if not hasattr(a, "sparsity"):
        return NotImplemented
    return type(a)(a)

def _np_diagflat(v, k=0):
    if not hasattr(v, "sparsity"):
        return NotImplemented
    # C-order flatten = column-major flatten of the transpose.
    flat = vec(v.T)
    n = flat.numel()
    if k == 0:
        return diag(flat)
    # Place `flat` on the k-th diagonal of an (n+|k|) x (n+|k|) matrix.
    size = n + abs(k)
    base = diag(flat)  # n x n
    if k > 0:
        # Shift right by k columns and pad bottom by k rows.
        top = horzcat(DM(n, k), base)            # n x size
        bot = DM(k, size)
        return vertcat(top, bot)
    else:
        # k < 0: shift down by |k| rows and pad right by |k| cols.
        bot = horzcat(base, DM(n, abs(k)))       # n x size
        top = DM(abs(k), size)
        return vertcat(top, bot)

def _np_empty_like(prototype, dtype=None, order='K', subok=True, shape=None):
    # dtype, order, subok are metadata-only for casadi's single value class.
    # shape= overrides the prototype shape: we can only honour it when the
    # prototype is non-symbolic (DM); for SX/MX we return NotImplemented and
    # let the fallback densify-and-numpy (which then returns an ndarray --
    # a reasonable answer for an explicit shape= request).
    if not hasattr(prototype, "sparsity"):
        return NotImplemented
    if shape is not None:
        # honour shape= by constructing a fresh dense zero matrix of the
        # given 2-D shape; reject 1-D / 3-D requests
        if isinstance(shape, int):
            r, c = shape, 1
        elif len(shape) == 1:
            r, c = int(shape[0]), 1
        elif len(shape) == 2:
            r, c = int(shape[0]), int(shape[1])
        else:
            return NotImplemented
        return type(prototype).zeros(r, c)
    return type(prototype).zeros(prototype.sparsity())

def _np_full(shape, fill_value, dtype=None, order='C', device=None,
             like=None):
    # dtype: ignored — casadi DM is float64.  For SX/MX, dtype is moot.
    # order: ignored — casadi has a single storage convention (CCS).
    # device: ignored — casadi has no device-placement concept.
    # like: ignored once dispatched.
    if device is not None or like is not None:
        return NotImplemented
    if not hasattr(fill_value, "sparsity"):
        return NotImplemented
    # numpy casadi values are always 2-D; we don't extend to 3-D+ shapes.
    if isinstance(shape, int):
        nr, nc = shape, 1
    else:
        s = tuple(shape)
        if len(s) == 1:
            nr, nc = int(s[0]), 1
        elif len(s) == 2:
            nr, nc = int(s[0]), int(s[1])
        else:
            return NotImplemented
    # fill_value must be a scalar; reject array-shaped values.
    if fill_value.numel() != 1:
        return NotImplemented
    return type(fill_value).ones(nr, nc) * fill_value

def _np_geomspace(start, stop, num=50, endpoint=True, dtype=None, axis=0):
    # dtype: metadata-only for casadi (no dtype concept) -- silently
    #        ignored only if None; otherwise hand off.
    # axis:  numpy places the new axis at `axis`. casadi values are
    #        always 2-D, so non-default axis -> NotImplemented.
    if dtype is not None:
        return NotImplemented
    if axis != 0:
        return NotImplemented
    n = int(num)
    if n <= 0:
        return NotImplemented
    # log-space then exponentiate. Equivalent to numpy's docstring:
    # "Similar to logspace, but with endpoints specified directly."
    log_start = log(start)
    log_stop = log(stop)
    if endpoint:
        lin = linspace(log_start, log_stop, n)
    else:
        full = linspace(log_start, log_stop, n + 1)
        lin = full[:n]
    return exp(lin)

def _np_logspace(start, stop, num=50, endpoint=True, base=10.0,
                 dtype=None, axis=0):
    # dtype: metadata-only; non-None -> hand off.
    # axis:  casadi values are 2-D so non-default axis -> NotImplemented.
    if dtype is not None:
        return NotImplemented
    if axis != 0:
        return NotImplemented
    n = int(num)
    if n <= 0:
        return NotImplemented
    if endpoint:
        lin = linspace(start, stop, n)
    else:
        full = linspace(start, stop, n + 1)
        lin = full[:n]
    # base ** lin; base may be a python scalar or a casadi value.
    return power(base, lin)

def _np_meshgrid(*xi, copy=True, sparse=False, indexing='xy'):
    # copy:   casadi values are reference-semantics; ignored.
    # sparse: numpy returns broadcastable 1-D-ish arrays. For our 2-D
    #         model this would still produce 2-D shapes, but semantics
    #         differ enough that we hand off when set.
    # indexing: 'xy' (default) vs 'ij'.
    if sparse:
        return NotImplemented
    if indexing not in ('xy', 'ij'):
        return NotImplemented
    if len(xi) != 2:
        # higher-D meshgrid not supported by 2-D casadi shape model.
        return NotImplemented
    x, y = xi
    # Flatten endpoints to 1-D casadi vectors (column).
    xv = vec(x)
    yv = vec(y)
    nx = xv.shape[0]
    ny = yv.shape[0]
    if indexing == 'xy':
        # numpy 'xy': X has shape (len(y), len(x)), Y has shape (len(y), len(x))
        X = repmat(xv.T, ny, 1)
        Y = repmat(yv, 1, nx)
    else:
        # 'ij': X shape (len(x), len(y)), Y shape (len(x), len(y))
        X = repmat(xv, 1, ny)
        Y = repmat(yv.T, nx, 1)
    return [X, Y]

def _np_tril(m, k=0):
    import numpy as np
    if not hasattr(m, "sparsity"):
        return NotImplemented
    nr, nc = m.shape
    # Build a dense mask: mask[i, j] = 1 if j - i <= k else 0
    # i.e. include entries on or below the k-th diagonal.
    if k == 0:
        return tril(m, True)
    if k == -1:
        return tril(m, False)
    # General k: explicit numeric mask.
    ii = np.arange(nr).reshape(-1, 1)
    jj = np.arange(nc).reshape(1, -1)
    mask = DM((jj - ii <= k).astype(float))
    return mask * m

def _np_triu(m, k=0):
    import numpy as np
    if not hasattr(m, "sparsity"):
        return NotImplemented
    nr, nc = m.shape
    if k == 0:
        return triu(m, True)
    if k == 1:
        return triu(m, False)
    # General k: mask[i,j] = 1 if j - i >= k else 0
    ii = np.arange(nr).reshape(-1, 1)
    jj = np.arange(nc).reshape(1, -1)
    mask = DM((jj - ii >= k).astype(float))
    return mask * m

def _np_vander(x, N=None, increasing=False):
    if not hasattr(x, "sparsity"):
        return NotImplemented
    # Treat x as a 1-D column vector.
    col = vec(x)
    M = col.numel()
    if N is None:
        N = M
    N = int(N)
    if N <= 0:
        # numpy returns shape (M, 0); casadi can't materialize 0-width.
        return NotImplemented
    powers = range(N) if increasing else range(N - 1, -1, -1)
    # Each column is col**p; stack horizontally.
    cols = [col ** p for p in powers]
    return horzcat(*cols)

def _np_isclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
    # equal_nan requires NaN comparison machinery; let numpy handle it
    # via the densify fallback (DM) or raise for symbolic.
    if equal_nan:
        return NotImplemented
    # Element-wise tolerance test.  fabs / le are CasADi free functions
    # so the result is a casadi expression of 0/1 values.
    return le(fabs(a - b), atol + rtol * fabs(b))

def _np_isrealobj(x):
    # CasADi values are always real-valued.
    return True

def _np_isfortran(a):
    # CasADi stores values in column-major (compressed-column-sparse)
    # layout -- F-order semantically.  Return True so the predicate
    # matches what the user gets out of casadi's vec(), reshape() etc.
    # (which are all column-major).  Without this entry, numpy's
    # np.isfortran tries to read `a.flags.fnc` and raises
    # AttributeError because DM/SX/MX don't expose .flags.
    return True

def _np_isscalar(element):
    # CasADi values are never Python scalars.  Numpy already returns
    # False for any object with __array_function__ etc., so this
    # handler is purely defensive.
    return False

def _np_around(a, decimals=0, out=None):
    # out= is not supported symbolically — route to fallback so DM
    # can densify and let numpy write into the user's buffer.
    if out is not None:
        return NotImplemented
    # casadi has no `round` free function; use sign+floor+fabs to
    # implement round-half-away-from-zero (tie-breaking differs from
    # numpy's banker's rounding at exact half-integers only).
    def _r(z):
        return _casadi.sign(z) * _casadi.floor(_casadi.fabs(z) + 0.5)
    if decimals == 0:
        return _r(a)
    # scale-round-unscale; works for any integer `decimals` (positive or
    # negative), e.g. decimals=-2 rounds to the nearest 100.
    factor = 10.0 ** int(decimals)
    return _r(a * factor) / factor

def _np_ediff1d(ary, to_end=None, to_begin=None):
    # to_end:   value(s) appended to the result (after diffs)
    # to_begin: value(s) prepended to the result (before diffs)
    # Returns a column-vector CasADi value (always 2-D shape (k,1)).
    if not hasattr(ary, "shape"):
        return NotImplemented

    # Flatten to a column vector (numpy ediff1d always operates on flat).
    v = vec(ary)
    n = v.shape[0]
    if n < 1:
        return NotImplemented

    if n >= 2:
        d = v[1:n, 0] - v[0:n-1, 0]
    else:
        d = DM.zeros(0, 1)

    parts = []
    if to_begin is not None:
        tb = vec(to_begin) if hasattr(to_begin, "shape") else DM(to_begin)
        tb = vec(tb)
        parts.append(tb)
    parts.append(d)
    if to_end is not None:
        te = vec(to_end) if hasattr(to_end, "shape") else DM(to_end)
        te = vec(te)
        parts.append(te)

    if len(parts) == 1:
        return parts[0]
    return vcat(parts)

def _np_fix(x, out=None):
    # out= forces the result to be written into a numpy buffer; route
    # back to numpy via NotImplemented + fallback.
    if out is not None:
        return NotImplemented
    # fix(x) = sign(x) * floor(abs(x))
    # Equivalent and smoother for branch-free codegen:
    #   if x >= 0: floor(x)  else: ceil(x)
    return if_else(x >= 0, floor(x), ceil(x), False)

def _prod_axis(a, axis):
    # Use a balanced-tree fold (depth O(log n)) instead of a left-deep
    # chain -- same total multiplications but shallower symbolic graph,
    # which matters for SX/MX AD on large reductions.  Split via
    # horzsplit/vertsplit so the per-column / per-row slices come from
    # a single casadi op rather than n Python-level indexing ops.
    if axis == 1:
        n = a.shape[1]
        if n == 0:
            return DM.ones(a.shape[0], 1)
        return _pairwise_reduce(lambda u, v: u * v, _casadi.horzsplit(a, 1))
    m = a.shape[0]
    if m == 0:
        return DM.ones(1, a.shape[1])
    return _pairwise_reduce(lambda u, v: u * v, _casadi.vertsplit(a, 1))

def _np_nanprod(a, axis=None, dtype=None, out=None,
                keepdims=False, initial=None, where=True):
    # dtype/keepdims: metadata + casadi-always-2-D (ignored).
    # out / where: route through fallback when non-default.
    # initial: multiplicative starting value.
    if not hasattr(a, "shape"):
        return NotImplemented
    if out is not None:
        return NotImplemented
    if where is not True:
        return NotImplemented

    # Mask NaNs to 1 via if_else (x == x is False iff x is NaN).
    masked = if_else(a == a, a, 1)

    if axis is None:
        flat = vec(masked)
        n = flat.shape[0]
        if n == 0:
            r = DM(1.0)
        else:
            r = _pairwise_reduce(lambda u, v: u * v,
                                 _casadi.vertsplit(flat, 1))
    elif axis in (0, -2):
        r = _prod_axis(masked, 0)
    elif axis in (1, -1):
        r = _prod_axis(masked, 1)
    else:
        return NotImplemented

    if initial is not None:
        r = r * initial
    return r

def _np_nansum(a, axis=None, dtype=None, out=None,
               keepdims=False, initial=None, where=True):
    # dtype/keepdims accepted but ignored (metadata + casadi-always-2-D).
    # out: cannot support symbolically.
    # where: element mask; route through fallback when non-default.
    # initial: additive starting value
    if not hasattr(a, "shape"):
        return NotImplemented
    if out is not None:
        return NotImplemented
    if where is not True:
        return NotImplemented

    # Mask NaNs to 0 via if_else / fmin / fmax of (a == a)
    # In CasADi:  is_nan(x) <=> (x != x). The DM 'eq' op is element-wise.
    masked = if_else(a == a, a, 0)

    if axis is None:
        r = sum1(sum2(masked))
    elif axis in (0, -2):
        r = sum1(masked)
    elif axis in (1, -1):
        r = sum2(masked)
    else:
        return NotImplemented

    if initial is not None:
        r = r + initial
    return r

def _np_prod(a, axis=None, dtype=None, out=None,
             keepdims=False, initial=None, where=True):
    # axis=None: reduce over all entries (return scalar)
    # axis=int : reduce along that axis (0 or 1 for 2-D casadi)
    # dtype:   metadata-only; ignored for casadi
    # out:     would write into user buffer; we can't support symbolically
    # keepdims: casadi values are always 2-D so the resulting shape
    #          already matches keepdims=True; for keepdims=False with
    #          axis=None numpy returns a 0-D scalar - we return (1,1).
    # initial: multiplicative starting value
    # where:   element mask (use 1 instead of element where mask is False)
    if not hasattr(a, "shape"):
        return NotImplemented
    if out is not None:
        return NotImplemented
    if where is not True:
        # symbolically expressible via if_else but not worth it here.
        return NotImplemented

    if axis is None:
        flat = vec(a)
        n = flat.shape[0]
        if n == 0:
            r = DM(1.0)
        else:
            r = _pairwise_reduce(lambda u, v: u * v,
                                 _casadi.vertsplit(flat, 1))
    elif axis in (0, -2):
        r = _prod_axis(a, 0)
    elif axis in (1, -1):
        r = _prod_axis(a, 1)
    else:
        return NotImplemented

    if initial is not None:
        r = r * initial
    return r

def _np_round(a, decimals=0, out=None):
    # out= forces the result to be written into a numpy buffer; route
    # back to numpy via NotImplemented + fallback.
    if out is not None:
        return NotImplemented
    if not isinstance(decimals, int):
        return NotImplemented
    # casadi has no `round` free function -- synthesize round-half-away-
    # from-zero (tie-breaking differs from numpy at exact half-integers).
    def _r(z):
        return _casadi.sign(z) * _casadi.floor(_casadi.fabs(z) + 0.5)
    if decimals == 0:
        return _r(a)
    # round to N decimals: scale, round, descale.
    scale = 10.0 ** decimals
    return _r(a * scale) / scale

def _np_sinc(x):
    # numpy.sinc has the signature numpy.sinc(x). No kwargs.
    px = pi * x
    # if_else handles the x=0 case symbolically; for AD, the derivative
    # at 0 is also defined (pi*cos(pi*x)/pi - sin(pi*x)/(pi*x^2)) but
    # if_else's branch will pick the 1.0 leaf there.
    return if_else(x == 0, 1.0, sin(px) / px)

def _to_const_list(test_elements):
    import numpy as np
    if isinstance(test_elements, (list, tuple)):
        flat = []
        for t in test_elements:
            sub = _to_const_list(t)
            if sub is None:
                return None
            flat.extend(sub)
        return flat
    if isinstance(test_elements, np.ndarray):
        return [float(v) for v in test_elements.ravel().tolist()]
    if hasattr(test_elements, "full"):
        # DM
        return [float(v) for v in test_elements.full().ravel().tolist()]
    if hasattr(test_elements, "is_constant") and test_elements.is_constant():
        return [float(v) for v in DM(test_elements).full().ravel().tolist()]
    # symbolic / unknown
    try:
        return [float(test_elements)]
    except Exception:
        return None

def _np_isin(element, test_elements, assume_unique=False, invert=False, kind=None):
    # assume_unique is a performance hint in numpy; no effect on result.
    # kind ('sort' / 'table' / None) is also a perf hint; ignore.
    if not hasattr(element, "shape"):
        return NotImplemented
    test_list = _to_const_list(test_elements)
    if test_list is None:
        return NotImplemented
    if len(test_list) == 0:
        # Empty test_elements: everything is (not) absent.
        zero = element * 0  # preserves type/shape
        return zero if not invert else zero + 1
    # OR of per-value equality masks via balanced-tree reduce
    # (depth O(log K) vs O(K) for a left-deep chain).
    acc = _pairwise_reduce(logic_or,
                           (eq(element, v) for v in test_list))
    if invert:
        acc = logic_not(acc)
    return acc


# --- batch 3: round-3 handlers ---

def _np_append(arr, values, axis=None):
    # axis=None: numpy flattens both inputs C-order (row-major), then
    #   concatenates.  CasADi has no 1-D and vec() is column-major, so
    #   transpose-then-vec to match numpy's element order.
    # axis=0:  vertical concat (rows must match column count).
    # axis=1:  horizontal concat.
    if axis is None:
        def _c_flat(x):
            if hasattr(x, "T"):
                return vec(x.T)
            return vec(DM(x).T)
        return vertcat(_c_flat(arr), _c_flat(values))
    if axis == 0:
        return vertcat(arr, values)
    if axis == 1:
        return horzcat(arr, values)
    return NotImplemented

def _np_array_split(ary, indices_or_sections, axis=0):
    if axis == 0:
        n = ary.shape[0]
        slicer = lambda lo, hi: ary[lo:hi, :]
    elif axis == 1:
        n = ary.shape[1]
        slicer = lambda lo, hi: ary[:, lo:hi]
    else:
        return NotImplemented

    if isinstance(indices_or_sections, int):
        # array_split: first (n%k) chunks get an extra row.
        k = indices_or_sections
        base, rem = divmod(n, k)
        sizes = [base + 1 if i < rem else base for i in range(k)]
        out = []
        cursor = 0
        for s in sizes:
            out.append(slicer(cursor, cursor + s))
            cursor += s
        return out
    # indices: explicit split points.
    try:
        cuts = [int(c) for c in indices_or_sections]
    except TypeError:
        return NotImplemented
    cuts = [0] + cuts + [n]
    return [slicer(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]

def _np_asfortranarray(a, dtype=None, like=None):
    import numpy as np
    # dtype: only float-like accepted; integer/bool casts would change
    #   semantics and CasADi can't represent them.
    # like: array-API protocol arg; only None makes sense.
    if like is not None:
        return NotImplemented
    if dtype is not None:
        try:
            np_dtype = np.dtype(dtype)
        except TypeError:
            return NotImplemented
        if not np.issubdtype(np_dtype, np.floating):
            return NotImplemented
    return a

def _np_block(arrays):
    # Single block: numpy returns it as-is (atleast_1d'd).
    if not isinstance(arrays, (list, tuple)):
        return arrays
    # Detect a flat row of blocks (no further nesting) -> hcat.
    flat_row = all(not isinstance(b, (list, tuple)) for b in arrays)
    if flat_row:
        return hcat(list(arrays))
    # Otherwise: outer list = rows; each row is a list of horizontal blocks.
    rows = []
    for r in arrays:
        if not isinstance(r, (list, tuple)):
            rows.append(r)
        else:
            # Disallow further nesting (would imply a 3rd axis).
            if any(isinstance(b, (list, tuple)) for b in r):
                return NotImplemented
            rows.append(hcat(list(r)))
    return vcat(rows)

def _np_broadcast_arrays(*args, subok=False):
    # subok: numpy-specific subclass passthrough; no analogue in CasADi.
    if subok:
        return NotImplemented
    # Determine target (row, col) via standard 2-D broadcasting.
    rows = 1
    cols = 1
    for a in args:
        sh = a.shape if hasattr(a, "shape") else (1, 1)
        r = sh[0] if len(sh) >= 1 else 1
        c = sh[1] if len(sh) >= 2 else 1
        if r != rows:
            if rows == 1:
                rows = r
            elif r != 1:
                return NotImplemented
        if c != cols:
            if cols == 1:
                cols = c
            elif c != 1:
                return NotImplemented
    out = []
    for a in args:
        sh = a.shape if hasattr(a, "shape") else (1, 1)
        r = sh[0] if len(sh) >= 1 else 1
        c = sh[1] if len(sh) >= 2 else 1
        out.append(repmat(a, rows // r, cols // c))
    return out

def _np_broadcast_to(array, shape, subok=False):
    # subok: numpy subclass passthrough; no analogue.
    if subok:
        return NotImplemented
    # Accept (m, n) tuple/list; we don't support 0-D or >2-D targets.
    sh = tuple(int(s) for s in shape)
    if len(sh) == 1:
        sh = (sh[0], 1)
    if len(sh) != 2:
        return NotImplemented
    m, n = sh
    a_sh = array.shape if hasattr(array, "shape") else (1, 1)
    r = a_sh[0] if len(a_sh) >= 1 else 1
    c = a_sh[1] if len(a_sh) >= 2 else 1
    # The target dim must be a multiple of source dim, with source dim
    # either 1 or equal.  (numpy is stricter; this matches the broadcast.)
    if r not in (1, m) or c not in (1, n):
        return NotImplemented
    return repmat(array, m // r, n // c)

def _np_delete(arr, obj, axis=None):
    if axis is None:
        # numpy: flatten arr, then delete by linear index. Round-tripping
        # through CasADi's column-flatten would change semantics relative
        # to numpy's row-flatten. Hand off.
        return NotImplemented
    if axis == 0:
        n = arr.shape[0]
        slicer = lambda i: arr[i, :]
    elif axis == 1:
        n = arr.shape[1]
        slicer = lambda i: arr[:, i]
    else:
        return NotImplemented

    # Normalize `obj` to a list of integer indices to remove.
    if isinstance(obj, slice):
        rm = set(range(*obj.indices(n)))
    elif isinstance(obj, int):
        rm = {obj % n if obj < 0 else obj}
    else:
        try:
            rm = set(int(i) % n if int(i) < 0 else int(i) for i in obj)
        except TypeError:
            return NotImplemented
    keep = [i for i in range(n) if i not in rm]
    if not keep:
        # Empty result; CasADi can represent (0, n) or (n, 0).
        if axis == 0:
            return arr[0:0, :]
        return arr[:, 0:0]
    parts = [slicer(i) for i in keep]
    return vcat(parts) if axis == 0 else hcat(parts)

def _np_flip(m, axis=None):
    if not hasattr(m, "shape"):
        return NotImplemented
    if axis is None:
        return m[::-1, ::-1]
    if isinstance(axis, (tuple, list)):
        axes = tuple(int(a) % 2 for a in axis)
    else:
        axes = (int(axis) % 2,)
    row = slice(None, None, -1) if 0 in axes else slice(None)
    col = slice(None, None, -1) if 1 in axes else slice(None)
    return m[row, col]

def _np_fliplr(m):
    if not hasattr(m, "shape"):
        return NotImplemented
    return m[:, ::-1]

def _np_flipud(m):
    if not hasattr(m, "shape"):
        return NotImplemented
    return m[::-1, :]

def _np_insert(arr, obj, values, axis=None):
    import numpy as np
    # Only support 2-D casadi with explicit axis = 0 or 1, plus integer
    # or list-of-integer obj.  Other cases fall back.
    if not hasattr(arr, "shape"):
        return NotImplemented
    if axis is None:
        # flatten-then-insert: would change shape rank, not 2-D-friendly.
        return NotImplemented
    if isinstance(obj, slice):
        return NotImplemented
    # Normalize obj to a sorted list of insertion positions.
    if isinstance(obj, (int, np.integer)):
        positions = [int(obj)]
    elif isinstance(obj, (list, tuple)):
        try:
            positions = [int(v) for v in obj]
        except Exception:
            return NotImplemented
    elif isinstance(obj, np.ndarray):
        try:
            positions = [int(v) for v in obj.ravel().tolist()]
        except Exception:
            return NotImplemented
    else:
        return NotImplemented

    nrows, ncols = arr.shape
    if axis == 0:
        n_along = nrows
    elif axis == 1:
        n_along = ncols
    else:
        return NotImplemented

    # numpy allows obj == n_along (append).  Normalize negative indices.
    positions = [p + n_along if p < 0 else p for p in positions]
    if any(p < 0 or p > n_along for p in positions):
        return NotImplemented

    # Build a permutation of slices: collect, in order, the original
    # rows/cols and the inserted rows/cols.
    n_inserts = len(positions)
    # Broadcast `values` against the missing-axis length.
    if axis == 0:
        # each insert is a row of length ncols
        vs = values
        if hasattr(vs, "shape"):
            # promote to 2-D row(s)
            pass
        # Build list of pieces along axis 0
        pieces = []
        last = 0
        # positions are not required to be sorted; numpy inserts AFTER
        # sorting them with stable order (each "before pos[i]" in the
        # ORIGINAL array).
        order = sorted(range(n_inserts), key=lambda i: positions[i])
        # Build a row-insertion plan in original-index order.
        plan = []  # list of (orig_pos, value_index_or_None)
        for k, i in enumerate(order):
            plan.append((positions[i], i))
        # Sweep original 0..nrows, dropping inserts at the right spot.
        ptr = 0
        for orig_pos, val_idx in plan:
            if orig_pos > last:
                pieces.append(arr[last:orig_pos, :])
                last = orig_pos
            # add one inserted row
            if hasattr(vs, "shape") and vs.shape == (n_inserts, ncols):
                pieces.append(vs[val_idx, :])
            else:
                # scalar or 1-D broadcast
                pieces.append(DM(vs) * DM.ones(1, ncols) if not hasattr(vs, "shape") else reshape(vs, 1, ncols))
        if last < nrows:
            pieces.append(arr[last:, :])
        return vcat(pieces)
    else:  # axis == 1
        pieces = []
        last = 0
        order = sorted(range(n_inserts), key=lambda i: positions[i])
        plan = [(positions[i], i) for i in order]
        vs = values
        for orig_pos, val_idx in plan:
            if orig_pos > last:
                pieces.append(arr[:, last:orig_pos])
                last = orig_pos
            if hasattr(vs, "shape") and vs.shape == (nrows, n_inserts):
                pieces.append(vs[:, val_idx])
            else:
                pieces.append(DM(vs) * DM.ones(nrows, 1) if not hasattr(vs, "shape") else reshape(vs, nrows, 1))
        if last < ncols:
            pieces.append(arr[:, last:])
        return hcat(pieces)

def _np_moveaxis(a, source, destination):
    import numpy as np
    if not hasattr(a, "shape"):
        return NotImplemented
    if isinstance(source, (tuple, list)):
        src = tuple(int(s) % 2 for s in source)
    else:
        src = (int(source) % 2,)
    if isinstance(destination, (tuple, list)):
        dst = tuple(int(d) % 2 for d in destination)
    else:
        dst = (int(destination) % 2,)
    if len(src) != len(dst):
        return NotImplemented
    # Build the resulting permutation of [0,1].
    perm = [0, 1]
    # Start with [0,1], move src[i] to dst[i].
    # Easiest: compute final position via numpy semantics on a probe.
    probe = np.array([0, 1])
    try:
        out = np.moveaxis(probe.reshape(2, 1), src, dst)
        # We just check which "slot" axis 0 ended up in.
    except Exception:
        return NotImplemented
    # Compute the permutation directly: which input axis is at each output axis.
    # Use numpy on a tagged shape.
    tag = np.empty((2, 3))  # shape distinguishes axes
    try:
        moved = np.moveaxis(tag, src, dst)
    except Exception:
        return NotImplemented
    # moved.shape tells us the permutation: output axis 0 has size
    # tag.shape[perm0], output axis 1 has size tag.shape[perm1].
    shape_in = (2, 3)
    perm = []
    for s in moved.shape:
        perm.append(shape_in.index(s))
        # mark used by replacing with sentinel
        shape_in = tuple(-1 if i == perm[-1] else v for i, v in enumerate(shape_in))
    if perm == [0, 1]:
        return a + 0
    if perm == [1, 0]:
        return a.T
    return NotImplemented

def _np_ndim(a):
    if not hasattr(a, "shape"):
        return NotImplemented
    return 2

def _np_ravel(a, order='C'):
    if not hasattr(a, "shape"):
        return NotImplemented
    o = str(order).upper()
    if o == 'F':
        return vec(a)
    if o in ('C', 'A', 'K'):
        return vec(a.T)
    return NotImplemented

def _np_rollaxis(a, axis, start=0):
    import numpy as np
    if not hasattr(a, "shape"):
        return NotImplemented
    # Use numpy on a tagged shape to compute the resulting permutation.
    try:
        moved = np.rollaxis(np.empty((2, 3)), int(axis), int(start))
    except Exception:
        return NotImplemented
    if moved.shape == (2, 3):
        return a + 0
    if moved.shape == (3, 2):
        return a.T
    return NotImplemented

def _np_rot90(m, k=1, axes=(0, 1)):
    if not hasattr(m, "shape"):
        return NotImplemented
    ax = tuple(int(a) for a in axes)
    if ax == (0, 1):
        sign = 1
    elif ax == (1, 0):
        sign = -1
    else:
        return NotImplemented
    kk = (sign * int(k)) % 4
    if kk == 0:
        return m + 0  # return a copy-like expression
    if kk == 1:
        return m.T[::-1, :]
    if kk == 2:
        return m[::-1, ::-1]
    # kk == 3
    return m[::-1, :].T

def _np_shape(a):
    if not hasattr(a, "shape"):
        return NotImplemented
    s = a.shape
    return (int(s[0]), int(s[1]))

def _np_size(a, axis=None):
    if not hasattr(a, "shape"):
        return NotImplemented
    nrow, ncol = int(a.shape[0]), int(a.shape[1])
    if axis is None:
        return nrow * ncol
    ax = int(axis) % 2
    return nrow if ax == 0 else ncol

def _np_squeeze(a, axis=None):
    if not hasattr(a, "shape"):
        return NotImplemented
    nrow, ncol = a.shape
    if axis is None:
        return a
    if isinstance(axis, (tuple, list)):
        axes = tuple(int(x) for x in axis)
    else:
        axes = (int(axis),)
    for ax in axes:
        ax = ax % 2  # normalize negatives
        if ax == 0 and nrow != 1:
            return NotImplemented
        if ax == 1 and ncol != 1:
            return NotImplemented
    return a

def _np_swapaxes(a, axis1, axis2):
    if not hasattr(a, "shape"):
        return NotImplemented
    ax1, ax2 = int(axis1) % 2, int(axis2) % 2
    if (ax1, ax2) in ((0, 0), (1, 1)):
        return a + 0
    if (ax1, ax2) in ((0, 1), (1, 0)):
        return a.T
    return NotImplemented

def _np_diagonal(a, offset=0, axis1=0, axis2=1):
    # numpy.diagonal extracts the diagonal of a 2-D submatrix indexed by
    # (axis1, axis2). CasADi values are always 2-D, so axis1/axis2 must be
    # the canonical (0, 1) or (1, 0) pair.
    if (axis1, axis2) not in ((0, 1), (1, 0)):
        return NotImplemented
    if (axis1, axis2) == (1, 0):
        a = a.T
    if offset == 0:
        # casadi.diag(M) extracts the diagonal of a 2-D matrix as a column vec
        return diag(a)
    # Off-diagonal extraction: index into a 2-D matrix
    nrow, ncol = a.shape
    if offset > 0:
        n = min(nrow, ncol - offset)
        if n <= 0:
            return NotImplemented
        return vcat([a[i, i + offset] for i in range(n)])
    else:
        n = min(nrow + offset, ncol)
        if n <= 0:
            return NotImplemented
        return vcat([a[i - offset, i] for i in range(n)])

def _np_linalg_cond(x, p=None):
    import numpy as np
    # p=2, -2, None require SVD — not available symbolically in CasADi.
    if p in (None, 2, -2):
        return NotImplemented
    Ainv = inv(x)
    # Delegate to _np_linalg_norm so we get the correct INDUCED matrix-norm
    # semantics (e.g. p=1 -> max column sum, p=inf -> max row sum).  An
    # earlier version used casadi.norm_1 / casadi.norm_inf directly which
    # are entrywise vector norms over vec(A) and give wrong cond numbers
    # by a factor of `n`.  See linalg/linalg_cond/test_handler.py.
    nx = _np_linalg_norm(x, ord=p)
    ni = _np_linalg_norm(Ainv, ord=p)
    if nx is NotImplemented or ni is NotImplemented:
        return NotImplemented
    return nx * ni

def _np_linalg_lstsq(a, b, rcond=None):
    import numpy as np
    # rcond: only affects the rank/singular-value cutoffs which we don't
    # compute; ignored.
    #
    # numpy ALWAYS returns the 4-tuple, so we cannot just return the
    # solution.  The dispatch wrapper unpacks the tuple if the caller
    # does e.g. `x, _, _, _ = np.linalg.lstsq(A, b)`.  We fill the
    # placeholder slots with empty objects.
    x = mtimes(pinv(a), b)
    # Residuals, rank, singular values are not symbolically meaningful.
    # numpy returns:
    #   residuals: shape (K,) or () when underdetermined / rank < n
    #   rank: int
    #   s: 1-D ndarray of singular values
    # We return placeholders; if the caller indexes [0] they get x.
    return (x, DM([]), DM(min(a.shape)), DM.nan(min(a.shape), 1))

def _np_linalg_matrix_power(a, n):
    import numpy as np
    # n must be an integer (numpy raises TypeError otherwise).
    if not isinstance(n, (int, np.integer)):
        return NotImplemented
    n = int(n)
    rows = a.shape[0]
    if a.shape[0] != a.shape[1]:
        return NotImplemented
    if n == 0:
        return DM.eye(rows)
    if n < 0:
        a = inv(a)
        n = -n
    # Exponentiation by squaring
    result = None
    base = a
    while n > 0:
        if n & 1:
            result = base if result is None else mtimes(result, base)
        n >>= 1
        if n > 0:
            base = mtimes(base, base)
    return result

def _np_linalg_multi_dot(arrays, out=None):
    # out is buffer output — ignored (NEP-18 fallback will handle).
    if out is not None:
        return NotImplemented
    arrays = list(arrays)
    if len(arrays) == 0:
        return NotImplemented
    if len(arrays) == 1:
        return arrays[0]
    result = arrays[0]
    for a in arrays[1:]:
        result = mtimes(result, a)
    return result

def _np_linalg_pinv(a, rcond=None, hermitian=False, rtol=None):
    # rcond / rtol: cutoff for small singular values; CasADi's pinv
    # doesn't expose this; return NotImplemented when set non-default.
    # hermitian=: numpy uses an eigendecomposition when True; we ignore
    # since casadi's pinv is dense-LSQR-style regardless.
    if rcond is not None or rtol is not None:
        return NotImplemented
    return pinv(a)

def _np_linalg_qr(a, mode='reduced'):
    Q, R = qr(a)
    if mode == 'reduced':
        return (Q, R)
    if mode == 'r':
        return R
    if mode == 'complete':
        # Only equals 'reduced' for square matrices.  For tall/wide
        # matrices we'd need to extend Q to a full-rank basis; not
        # supported.
        m, n = a.shape
        if m == n:
            return (Q, R)
        return NotImplemented
    if mode == 'raw':
        return NotImplemented
    return NotImplemented

def _np_linalg_slogdet(a):
    # No kwargs in numpy.linalg.slogdet.
    d = det(a)
    return sign(d), log(fabs(d))

def _np_tensordot(a, b, axes=2):
    if isinstance(axes, int):
        if axes == 0:
            # outer product: a[i,j] * b[k,l] -- N-D result.  For
            # CasADi (always 2-D) we only do this if both are vectors,
            # giving the outer product matrix.
            if a.shape[1] == 1 and b.shape[1] == 1:
                return mtimes(a, b.T)
            return NotImplemented
        if axes == 1:
            # standard matmul of last axis of a with first axis of b.
            return mtimes(a, b)
        if axes == 2:
            # Frobenius inner: sum_{ij} a[i,j] * b[i,j]
            if a.shape == b.shape:
                return sum1(sum2(a * b))
            return NotImplemented
        return NotImplemented
    # axes = (axes_a, axes_b) tuple
    if isinstance(axes, (list, tuple)) and len(axes) == 2:
        aa, bb = axes
        if isinstance(aa, int):
            aa = [aa]
        if isinstance(bb, int):
            bb = [bb]
        if len(aa) != len(bb):
            return NotImplemented
        # 2-D only: axis 0 = rows, axis 1 = cols.
        if len(aa) == 1:
            ax_a, ax_b = aa[0], bb[0]
            ta = a if ax_a == 1 else a.T   # contracted axis to the right
            tb = b if ax_b == 0 else b.T   # contracted axis to the left
            return mtimes(ta, tb)
        if len(aa) == 2 and a.shape == b.shape:
            # full contraction
            if sorted(aa) == [0, 1] and sorted(bb) == [0, 1]:
                # may need transpose if axis order differs but for
                # equal-shape, result is the same scalar.
                return sum1(sum2(a * b))
        return NotImplemented
    return NotImplemented

def _np_vdot(a, b):
    # numpy.vdot flattens both arguments to 1-D and returns the dot product.
    # For complex it conjugates the first arg; CasADi is real-only so no conj.
    # CasADi values are always 2-D — vec() makes them column-stacked 1-D
    # (in column-major / Fortran order). numpy.vdot ravels in C order, so to
    # match numpy's ravel-to-1-D semantics for non-trivial 2-D inputs we
    # transpose first so vec() picks elements in the right order.
    av = vec(a.T) if hasattr(a, "T") else vec(a)
    bv = vec(b.T) if hasattr(b, "T") else vec(b)
    return dot(av, bv)

def _np_allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
    import numpy as np
    def _to_np(z):
        if hasattr(z, "is_constant"):
            if not z.is_constant():
                return None
            return z.full() if hasattr(z, "full") else np.asarray(DM(z).full())
        return np.asarray(z)
    a_np = _to_np(a)
    b_np = _to_np(b)
    if a_np is None or b_np is None:
        # symbolic non-constant: can't decide, fall back to NEP-18 default
        # (which will raise a TypeError on SX/MX -- the desired behaviour).
        return NotImplemented
    return bool(np.allclose(a_np, b_np, rtol=rtol, atol=atol, equal_nan=equal_nan))

def _np_array_equal(a1, a2, equal_nan=False):
    import numpy as np
    def _to_np(z):
        if hasattr(z, "is_constant"):
            if not z.is_constant():
                return None
            return z.full() if hasattr(z, "full") else np.asarray(DM(z).full())
        return np.asarray(z)
    a = _to_np(a1)
    b = _to_np(a2)
    if a is None or b is None:
        # symbolic non-constant; fallback densifies (errors on SX/MX).
        return NotImplemented
    return bool(np.array_equal(a, b, equal_nan=equal_nan))

def _np_iscomplexobj(x):
    # CasADi DM/SX/MX are always real-valued; complex dtype is not supported.
    return False

def _np_isneginf(x, out=None):
    import numpy as np
    # out=: numpy-only output array; not bridged.  When the user passes a
    # non-None out, return NotImplemented so the NEP-18 fallback densifies.
    if out is not None:
        return NotImplemented
    if hasattr(x, "is_constant") and not x.is_constant():
        return DM.zeros(x.sparsity())
    f = x.full() if hasattr(x, "full") else np.asarray(x)
    return DM(np.isneginf(f).astype(float))

def _np_isposinf(x, out=None):
    import numpy as np
    # out=: numpy-only output array; not bridged.  Pass-through to fallback.
    if out is not None:
        return NotImplemented
    if hasattr(x, "is_constant") and not x.is_constant():
        return DM.zeros(x.sparsity())
    f = x.full() if hasattr(x, "full") else np.asarray(x)
    return DM(np.isposinf(f).astype(float))

def _np_convolve(a, v, mode='full'):
    # CasADi values are 2-D. Flatten both to column vectors.
    a_flat = vec(a)
    v_flat = vec(v)
    N = a_flat.shape[0]
    M = v_flat.shape[0]
    L_full = N + M - 1

    # Result of full convolution at index k:
    #     y[k] = sum_{j=max(0,k-N+1)}^{min(M-1,k)} v[j] * a[k-j]
    # Pre-split both inputs into scalar pieces with one casadi op each
    # so the per-(k,j) accesses are plain python list indexing rather
    # than N*M repeated casadi indexing ops.  Inner sum uses balanced
    # tree reduce for O(log K) symbolic depth per output element.
    v_pieces = _casadi.vertsplit(v_flat, 1)
    a_pieces = _casadi.vertsplit(a_flat, 1)
    rows = []
    for k in range(L_full):
        j_lo = max(0, k - N + 1)
        j_hi = min(M - 1, k)
        terms = [v_pieces[j] * a_pieces[k - j] for j in range(j_lo, j_hi + 1)]
        if not terms:
            rows.append(0)
        else:
            rows.append(_pairwise_reduce(lambda u, w: u + w, terms))
    full = vcat(rows)

    if mode == 'full':
        return full
    if mode == 'same':
        out_len = max(N, M)
        start = (L_full - out_len) // 2
        return full[start:start + out_len]
    if mode == 'valid':
        out_len = max(N, M) - min(N, M) + 1
        start = min(N, M) - 1
        return full[start:start + out_len]
    raise ValueError("convolve: mode must be 'full', 'same', or 'valid'")

def _np_cumprod(a, axis=None, dtype=None, out=None):
    # dtype/out: numpy machinery; ignored — CasADi has fixed dtype and
    # no out-array semantics.
    if out is not None:
        return NotImplemented

    # Cumulative product is inherently chain-shaped (output[i] depends
    # on output[i-1]), so the chain depth O(n) is unavoidable here.
    # We can still cut N python-level indexings by splitting once.
    if axis is None:
        # numpy's axis=None flattens in C-order (row-major).
        # casadi.vec is column-major, so transpose first.
        flat = vec(a.T)
        if flat.shape[0] == 0:
            return DM.zeros(0, 1)
        pieces = _casadi.vertsplit(flat, 1)
        cols = [pieces[0]]
        for p in pieces[1:]:
            cols.append(cols[-1] * p)
        return vcat(cols)

    if axis == 0:
        # Cumulative product down each column.
        if a.shape[0] == 0:
            return DM.zeros(0, a.shape[1])
        pieces = _casadi.vertsplit(a, 1)
        rows = [pieces[0]]
        for p in pieces[1:]:
            rows.append(rows[-1] * p)
        return vcat(rows)

    if axis == 1:
        # Cumulative product across each row.
        if a.shape[1] == 0:
            return DM.zeros(a.shape[0], 0)
        pieces = _casadi.horzsplit(a, 1)
        cols = [pieces[0]]
        for p in pieces[1:]:
            cols.append(cols[-1] * p)
        return hcat(cols)

    return NotImplemented

def _to_const_list(v):
    import numpy as np
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    if isinstance(v, np.ndarray):
        return [float(x) for x in v.ravel().tolist()]
    if hasattr(v, "full"):
        return [float(x) for x in v.full().ravel().tolist()]
    if hasattr(v, "is_constant") and v.is_constant():
        return [float(x) for x in DM(v).full().ravel().tolist()]
    return None

def _np_interp(x, xp, fp, left=None, right=None, period=None):
    # Build a piecewise-linear chain by hand via if_else.  Why not use
    # casadi.interp1d?  Because its signature is
    #   interp1d([float] x, DM v, [float] xq, str mode, bool equidistant)
    # which requires xq to be a list of floats, ruling out symbolic
    # query points.  The hand-built chain supports SX/MX queries (and
    # numeric DM queries) uniformly.
    xp_list = _to_const_list(xp)
    fp_list = _to_const_list(fp)
    if xp_list is None or fp_list is None:
        return NotImplemented
    n = len(xp_list)
    if n != len(fp_list) or n < 2:
        return NotImplemented

    if period is not None:
        # Wrap x into [xp[0], xp[0]+period); after wrapping, numpy
        # ignores left/right per the documented spec.
        x0 = xp_list[0]
        x = x0 + fmod(fmod(x - x0, period) + period, period)
        left = None
        right = None

    # numpy clamps to fp[0] / fp[-1] outside the range unless left/right given.
    lv = fp_list[0] if left is None else left
    rv = fp_list[-1] if right is None else right

    # Build from the right: y starts as the "above-the-range" value.
    y = rv
    for i in range(n - 2, -1, -1):
        x0, x1 = xp_list[i], xp_list[i + 1]
        y0, y1 = fp_list[i], fp_list[i + 1]
        # Linear segment between (x0,y0) and (x1,y1).  Constants cancel
        # so this is a single mul+add on the symbolic x.
        seg = y0 + (x - x0) * ((y1 - y0) / (x1 - x0))
        y = if_else(x < x1, seg, y)
    # Below-range clamp / left override.
    y = if_else(x < xp_list[0], lv, y)
    return y

def _np_nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None):
    import numpy as np
    # `copy` is irrelevant: CasADi expressions are immutable (always a
    # fresh node) so the in-place fast path doesn't apply.
    if not hasattr(x, "shape"):
        return NotImplemented

    # Match numpy's defaults for posinf / neginf (largest/smallest
    # representable finite double).
    if posinf is None:
        posinf = np.finfo(np.float64).max
    if neginf is None:
        neginf = np.finfo(np.float64).min

    INF = float('inf')
    # Order matters: do NaN first (via self-equality), then infinities.
    y = if_else(x == x, x, nan)
    y = if_else(y == INF, posinf, y)
    y = if_else(y == -INF, neginf, y)
    return y

def _np_nancumsum(a, axis=None, dtype=None, out=None):
    # dtype: metadata, ignored.
    # out:   cannot support symbolically.
    if out is not None:
        return NotImplemented
    if not hasattr(a, "shape"):
        return NotImplemented

    # Mask NaNs to 0 via if_else / self-equality (NaN != NaN).
    masked = if_else(a == a, a, 0)

    # Cumulative sum: casadi.cumsum(x, axis) exists.
    if axis is None:
        return cumsum(vec(masked), 0)
    return cumsum(masked, int(axis))

def _np_nanmax(a, axis=None, out=None, keepdims=False, initial=None, where=True):
    # axis=None: reduce over all elements
    # axis=0/1: row/column reduction
    # axis=tuple: not natively supported by mmin/mmax; NotImplemented
    # out=: ignored / unsupported for CasADi
    # keepdims=: CasADi values are always 2-D; default True-like behaviour
    # initial=: prepended into the reduction; supported scalar
    # where=: would mask elements; not supported symbolically
    if out is not None:
        return NotImplemented
    if where is not True:
        return NotImplemented
    if axis is None:
        r = mmax(a)
        if initial is not None:
            r = fmax(r, initial)
        return r
    if isinstance(axis, tuple):
        return NotImplemented
    if axis == 0:
        # Column-wise maxima -> (1, ncols).  Balanced-tree fmax fold
        # has O(log n_rows) symbolic depth vs O(n_rows) for a chain.
        # Use vertsplit so the per-row slices come from one casadi op
        # rather than n_rows Python-level a[i,:] indexings.
        r = _pairwise_reduce(fmax, _casadi.vertsplit(a, 1))
        if initial is not None:
            r = fmax(r, initial)
        return r
    if axis == 1:
        # Row-wise maxima -> (nrows, 1).
        r = _pairwise_reduce(fmax, _casadi.horzsplit(a, 1))
        if initial is not None:
            r = fmax(r, initial)
        return r
    return NotImplemented

def _np_nanmin(a, axis=None, out=None, keepdims=False, initial=None, where=True):
    # keepdims: ignored (casadi always-2-D).
    # out:      cannot support symbolically.
    # initial:  upper bound on result; mappable but rare.
    # where:    element mask; route through fallback when non-default.
    if not hasattr(a, "shape"):
        return NotImplemented
    if out is not None or where is not True or initial is not None:
        return NotImplemented
    if axis is not None:
        return NotImplemented

    # Mask NaNs to +inf so they don't influence the min.
    INF = float('inf')
    masked = if_else(a == a, a, INF)
    return mmin(masked)

def _np_real_if_close(a, tol=100):
    # tol: meaningless for real-only CasADi; accept and ignore.
    return a + 0  # cheap copy preserves type

def _np_average(a, axis=None, weights=None, returned=False, keepdims=False):
    # returned: would return (avg, sum_of_weights) tuple — non-default falls back.
    # keepdims: casadi is always 2-D; reduction along an axis already keeps the other.
    if returned:
        return NotImplemented
    n_rows, n_cols = a.shape
    if weights is None:
        if axis is None:
            n = n_rows * n_cols
            return sum2(sum1(a)) / n
        if axis == 0:
            return sum1(a) / n_rows
        if axis == 1:
            return sum2(a) / n_cols
        return NotImplemented
    # Weighted average: numpy requires weights.shape == a.shape (along axis).
    w = weights
    if axis is None:
        return sum2(sum1(a * w)) / sum2(sum1(w))
    if axis == 0:
        # weights is 1-D of length n_rows in numpy; promote to (n_rows,1) and broadcast.
        if hasattr(w, "shape") and len(getattr(w, "shape", ())) == 1:
            w = reshape(w, n_rows, 1)
            w = repmat(w, 1, n_cols)
        return sum1(a * w) / sum1(w)
    if axis == 1:
        if hasattr(w, "shape") and len(getattr(w, "shape", ())) == 1:
            w = reshape(w, 1, n_cols)
            w = repmat(w, n_rows, 1)
        return sum2(a * w) / sum2(w)
    return NotImplemented

def _np_corrcoef(x, y=None, rowvar=True, bias=None, ddof=None, dtype=None):
    # dtype, bias, ddof: bias/ddof are deprecated no-ops in numpy >= 1.10.
    # y: stacking-with-y branch — fall back to numpy.
    if y is not None:
        return NotImplemented
    X = x if rowvar else x.T
    nv, nobs = X.shape
    mu = sum2(X) / nobs                       # (nv, 1)
    Xc = X - repmat(mu, 1, nobs)
    C = mtimes(Xc, Xc.T) / (nobs - 1)         # ddof choice cancels in normalization
    d = sqrt(diag(C))                          # (nv, 1)
    norm = mtimes(d, d.T)
    return C / norm

def _np_cov(m, y=None, rowvar=True, bias=False, ddof=None, fweights=None,
            aweights=None, dtype=None):
    # dtype: metadata, ignored.
    # y, fweights, aweights: non-trivial, fall back to numpy.
    if y is not None or fweights is not None or aweights is not None:
        return NotImplemented
    # Orient so rows = variables, cols = observations.
    X = m if rowvar else m.T
    nv, nobs = X.shape
    # ddof default: 1 if bias else 0  (using numpy's bias→ddof=0, unbiased→ddof=1).
    if ddof is None:
        ddof = 0 if bias else 1
    # Center.
    mu = sum2(X) / nobs                 # (nv, 1)
    Xc = X - repmat(mu, 1, nobs)
    return mtimes(Xc, Xc.T) / (nobs - ddof)

def _np_digitize(x, bins, right=False):
    # bins must be constant: the partition is structural, not symbolic.
    edges = _to_const_list(bins)
    if edges is None:
        return NotImplemented
    n = len(edges)
    if n == 0:
        return x * 0  # everything in "bin 0"
    # Detect monotonicity (numpy supports both ascending and descending).
    ascending = all(edges[i] <= edges[i+1] for i in range(n-1))
    descending = all(edges[i] >= edges[i+1] for i in range(n-1))
    if not (ascending or descending):
        return NotImplemented
    # For ascending bins, result = sum_i (cond_i), where cond_i is:
    #   right=False:  x >= edges[i]
    #   right=True :  x >  edges[i]
    # For descending bins, returns i such that bins[i-1] (>|>=) x (>=|>) bins[i].
    # The count-of-passed-edges trick swaps the inequality DIRECTION ("past"
    # means smaller value) AND swaps which kind of tie counts as "passed",
    # i.e. descending-right=False uses strict-less-than while
    # descending-right=True uses non-strict less-than-or-equal.  An earlier
    # version had the descending pair flipped, miscounting ties.
    acc = x * 0  # zero with proper type/shape
    for e in edges:
        if ascending:
            cond = (x > e) if right else (x >= e)
        else:
            cond = (x <= e) if right else (x < e)
        acc = acc + cond
    return acc

def _np_mean(a, axis=None, dtype=None, out=None, keepdims=False, where=True):
    # dtype, out: metadata / write-target — no-op for casadi (always dense, no dtype).
    # keepdims: casadi is always 2-D, so axis-reductions already keep the
    #   "other" dim (sum1 -> 1xN, sum2 -> Nx1); we don't need an extra branch.
    # where: data-dependent boolean mask. Defer to numpy fallback if non-default.
    if where is not True:
        return NotImplemented
    n_rows, n_cols = a.shape
    if axis is None:
        n = n_rows * n_cols
        return sum2(sum1(a)) / n
    if axis == 0:
        return sum1(a) / n_rows
    if axis == 1:
        return sum2(a) / n_cols
    return NotImplemented

def _np_ptp(a, axis=None, out=None, keepdims=None):
    # out, keepdims: see _np_mean rationale.
    # casadi has no per-axis mmax/mmin; only axis=None is supported.
    if axis is None:
        return mmax(a) - mmin(a)
    return NotImplemented

def _np_std(a, axis=None, dtype=None, out=None, ddof=0, keepdims=False,
            where=True, mean=None, correction=None):
    if where is not True or mean is not None or correction is not None:
        return NotImplemented
    n_rows, n_cols = a.shape
    if axis is None:
        n = n_rows * n_cols
        mu = sum2(sum1(a)) / n
        d = a - mu
        return sqrt(sum2(sum1(d * d)) / (n - ddof))
    if axis == 0:
        mu = sum1(a) / n_rows
        d = a - repmat(mu, n_rows, 1)
        return sqrt(sum1(d * d) / (n_rows - ddof))
    if axis == 1:
        mu = sum2(a) / n_cols
        d = a - repmat(mu, 1, n_cols)
        return sqrt(sum2(d * d) / (n_cols - ddof))
    return NotImplemented

def _np_var(a, axis=None, dtype=None, out=None, ddof=0, keepdims=False,
            where=True, mean=None, correction=None):
    # dtype, out, keepdims: see _np_mean rationale.
    # where, mean (precomputed mean override), correction (alias for ddof):
    #   non-trivial branches — fall back to numpy when set.
    if where is not True or mean is not None or correction is not None:
        return NotImplemented
    n_rows, n_cols = a.shape
    if axis is None:
        n = n_rows * n_cols
        mu = sum2(sum1(a)) / n
        d = a - mu
        return sum2(sum1(d * d)) / (n - ddof)
    if axis == 0:
        mu = sum1(a) / n_rows           # 1xN
        d = a - repmat(mu, n_rows, 1)
        return sum1(d * d) / (n_rows - ddof)
    if axis == 1:
        mu = sum2(a) / n_cols           # Nx1
        d = a - repmat(mu, 1, n_cols)
        return sum2(d * d) / (n_cols - ddof)
    return NotImplemented


# --- batch 4: round-4 handlers ---

def _expand_pad_width(pad_width, ndim=2):
    # Returns a list of (before, after) per axis, length == ndim.
    if isinstance(pad_width, int):
        return [(pad_width, pad_width)] * ndim
    pw = list(pad_width)
    if len(pw) == 2 and all(isinstance(x, int) for x in pw):
        # (before, after) shared across all axes
        return [(pw[0], pw[1])] * ndim
    if len(pw) == ndim:
        out = []
        for entry in pw:
            entry = list(entry)
            if len(entry) == 1:
                out.append((entry[0], entry[0]))
            elif len(entry) == 2:
                out.append((entry[0], entry[1]))
            else:
                return None
        return out
    return None

def _expand_const_values(constant_values, ndim=2):
    # Returns a list of (before_val, after_val) per axis, length == ndim.
    if isinstance(constant_values, (int, float)):
        return [(constant_values, constant_values)] * ndim
    cv = list(constant_values)
    if len(cv) == 2 and all(isinstance(x, (int, float)) for x in cv):
        return [(cv[0], cv[1])] * ndim
    if len(cv) == ndim:
        out = []
        for entry in cv:
            entry = list(entry)
            if len(entry) == 1:
                out.append((entry[0], entry[0]))
            elif len(entry) == 2:
                out.append((entry[0], entry[1]))
            else:
                return None
        return out
    return None

def _np_pad(array, pad_width, mode='constant', constant_values=0):
    # Only 'constant' mode is supported symbolically.  Other modes
    # ('edge', 'reflect', 'linear_ramp', 'maximum', 'mean', 'median',
    # 'minimum', 'symmetric', 'wrap', 'empty', or a callable) need
    # value-dependent or wrap-around access patterns that don't have
    # clean symbolic equivalents — return NotImplemented and let
    # the dispatcher fall back (which works for DM and errors clearly
    # for SX/MX).
    if mode != 'constant':
        return NotImplemented
    if not hasattr(array, "shape"):
        return NotImplemented
    pw = _expand_pad_width(pad_width, ndim=2)
    if pw is None:
        return NotImplemented
    cv = _expand_const_values(constant_values, ndim=2)
    if cv is None:
        return NotImplemented
    (b0, a0), (b1, a1) = pw
    (cb0, ca0), (cb1, ca1) = cv
    # numpy.pad pads axis-by-axis in order 0, 1, ..., n-1.  That order
    # matters when each axis has its own constant: the new rows added by
    # the axis-0 step get extended by axis-1 padding, and the axis-1
    # constants fill the resulting corner cells (rather than the axis-0
    # constants -- which would happen if we padded axis 1 first).
    cls = type(array)
    cols = array.shape[1]
    pieces0 = []
    if b0 > 0:
        pieces0.append(cls(b0, cols) + cb0)
    pieces0.append(array)
    if a0 > 0:
        pieces0.append(cls(a0, cols) + ca0)
    mid = vertcat(*pieces0) if len(pieces0) > 1 else pieces0[0]
    rows = mid.shape[0]
    pieces1 = []
    if b1 > 0:
        pieces1.append(cls(rows, b1) + cb1)
    pieces1.append(mid)
    if a1 > 0:
        pieces1.append(cls(rows, a1) + ca1)
    return horzcat(*pieces1) if len(pieces1) > 1 else pieces1[0]

def _np_resize(a, new_shape):
    if not hasattr(a, "shape"):
        return NotImplemented
    # new_shape can be an int or a tuple/list.
    if isinstance(new_shape, int):
        # 1-D request: not representable in casadi's 2-D model.
        # Use (1, new_shape) as a natural 2-D promotion.
        nr, nc = 1, new_shape
    else:
        ns = list(new_shape)
        if len(ns) == 1:
            nr, nc = 1, int(ns[0])
        elif len(ns) == 2:
            nr, nc = int(ns[0]), int(ns[1])
        else:
            # higher-rank shapes don't fit casadi's 2-D model
            return NotImplemented
    old_total = int(a.shape[0]) * int(a.shape[1])
    new_total = nr * nc
    if old_total == 0:
        return NotImplemented
    # Row-major flatten == column-major flatten of the transpose.
    flat = vec(a.T)              # (old_total, 1) column-vector, row-major order
    # Tile enough times to cover new_total.
    n_tiles = (new_total + old_total - 1) // old_total
    tiled = repmat(flat, n_tiles, 1)
    sliced = tiled[:new_total, 0]
    # Reshape to (nc, nr) column-major == (nr, nc) row-major.
    return reshape(sliced, nc, nr).T

def _np_split(ary, indices_or_sections, axis=0):
    # axis must be 0 or 1; anything else is meaningless for 2-D casadi.
    if axis not in (0, 1):
        return NotImplemented
    if not hasattr(ary, "shape"):
        return NotImplemented
    n = ary.shape[axis]
    # Build offsets using the same helper used by _np_hsplit / _np_vsplit.
    offsets = _np_split_offsets(n, indices_or_sections)
    if offsets is None:
        return NotImplemented
    return vertsplit(ary, offsets) if axis == 0 else horzsplit(ary, offsets)

def _np_emath_arccos(x):
    # No kwargs in emath.arccos.
    return arccos(x)

def _np_emath_arcsin(x):
    # emath.arcsin takes no kwargs.
    return arcsin(x)

def _np_emath_arctanh(x):
    # No kwargs.
    return arctanh(x)

def _np_emath_log2(x):
    # No kwargs.
    return log(x) / log(2)

def _np_emath_logn(n, x):
    # No kwargs.
    return log(x) / log(n)

def _np_emath_power(x, p):
    # No kwargs.
    return power(x, p)

def _np_choose(a, choices, out=None, mode='raise'):
    # out unsupported (immutable casadi values).
    if out is not None:
        return NotImplemented
    # mode 'raise' is the only one whose semantics match a direct lookup.
    if mode != 'raise':
        return NotImplemented
    n = len(choices)
    if n == 0:
        return NotImplemented
    # Fold from last index down to first: a == 0 ? choices[0] :
    #   (a == 1 ? choices[1] : ...)
    result = choices[n - 1]
    for i in range(n - 2, -1, -1):
        result = if_else(a == i, choices[i], result)
    return result

def _np_fill_diagonal(a, val, wrap=False):
    # wrap=True (block-wrap for tall matrices) — not implemented; fall back.
    if wrap:
        return NotImplemented
    n, m = a.shape
    k = min(n, m)
    # Diagonal index list:
    idx = list(range(k))
    # Build a copy with diagonal replaced.  val may be scalar or length-k.
    val_arr = val
    try:
        val_arr = DM(val)
    except Exception:
        pass
    out = a + 0  # force a fresh expression
    if hasattr(val_arr, 'is_scalar') and val_arr.is_scalar():
        for i in idx:
            out[i, i] = val_arr
    else:
        for j, i in enumerate(idx):
            out[i, i] = val_arr[j]
    # numpy returns None (in-place).  Return the new matrix instead —
    # users wanting the numpy semantics should use fallback.
    return out

def _np_select(condlist, choicelist, default=0):
    if len(condlist) != len(choicelist):
        return NotImplemented
    if len(condlist) == 0:
        return NotImplemented
    # Build a right-to-left chain: if_else(c0, ch0, if_else(c1, ch1, ... default))
    result = default
    for cond, choice in zip(reversed(list(condlist)), reversed(list(choicelist))):
        result = if_else(cond, choice, result)
    return result

def _np_take(a, indices, axis=None, out=None, mode='raise'):
    import numpy as np
    # out: only accept the default; otherwise hand back to numpy fallback.
    if out is not None:
        return NotImplemented
    # mode: 'raise' is the only one whose semantics match a straight slice;
    # 'wrap' / 'clip' would change which element is fetched.
    if mode != 'raise':
        return NotImplemented
    # indices must be a concrete integer sequence — we can't index by a
    # casadi symbol along an arbitrary axis.
    try:
        idx = np.asarray(indices, dtype=int)
    except Exception:
        return NotImplemented

    if axis is None:
        # Per numpy semantics: a is flattened (C-order) before indexing.
        # CasADi storage is column-major; emulate C-order flatten via
        # row-major reshape using transpose.
        flat = vec(a.T)   # row-major flatten
        flat_idx = idx.ravel().tolist()
        out_val = flat[flat_idx, 0]
        # Reshape back to indices' shape (1-D result if idx is 1-D).
        if idx.ndim <= 1:
            return out_val
        # 2-D indices → reshape to (rows, cols) in row-major sense.
        rows, cols = idx.shape
        return reshape(out_val, cols, rows).T

    if axis == 0:
        return a[idx.ravel().tolist(), :]
    if axis == 1:
        return a[:, idx.ravel().tolist()]
    # Higher axes don't apply to 2-D casadi values.
    return NotImplemented

def _np_tril_indices_from(arr, k=0):
    import numpy as np
    # arr: CasADi value. shape is static.
    # k: diagonal offset.
    n, m = arr.shape
    return np.tril_indices(n, k=k, m=m)

def _np_triu_indices_from(arr, k=0):
    import numpy as np
    # arr supplies shape only.
    n, m = arr.shape
    if n != m:
        # numpy actually only supports 2-D arrays of any shape; the
        # constraint is just ndim==2.
        pass
    return np.triu_indices(n, k=k, m=m)

def _np_array_equiv(a1, a2):
    import numpy as np
    def _to_np(z):
        if hasattr(z, "is_constant"):
            if not z.is_constant():
                return None
            return z.full() if hasattr(z, "full") else np.asarray(DM(z).full())
        return np.asarray(z)
    a = _to_np(a1)
    b = _to_np(a2)
    if a is None or b is None:
        # symbolic non-constant: cannot decide; fall back (errors on SX/MX).
        return NotImplemented
    return bool(np.array_equiv(a, b))

def _np_iscomplex(x):
    # CasADi DM/SX/MX are always real-valued; every element is non-complex.
    # Return a zero-valued same-shape expression rather than a python False,
    # to match numpy.iscomplex's element-wise return-shape contract.
    if hasattr(x, "sparsity"):
        return type(x).zeros(x.sparsity().size1(), x.sparsity().size2())
    return NotImplemented

def _np_isreal(x):
    # CasADi DM/SX/MX are always real-valued; every element is real.
    if hasattr(x, "sparsity"):
        return type(x).ones(x.sparsity().size1(), x.sparsity().size2())
    return NotImplemented

def _np_count_nonzero(a, axis=None, keepdims=False):
    # axis=None  -> flatten and reduce to scalar
    # axis=0     -> sum1 (column-wise, shape (1, ncol))
    # axis=1     -> sum2 (row-wise,    shape (nrow, 1))
    # keepdims=False: axis=None gives scalar; axis=0/1 gives a 2-D
    #   shape with a singleton; in CasADi everything is already 2-D
    #   so the (1, n) / (n, 1) outputs are already the "keepdims" form.
    #   For keepdims=False on axis=0 numpy returns shape (n,) which is
    #   incompatible with CasADi's 2-D model; we keep the 2-D shape
    #   either way (acceptable: most consumers re-broadcast).

    mask = (a != 0)

    if axis is None:
        return sum1(sum2(mask))

    if axis == 0:
        return sum1(mask)
    if axis == 1:
        return sum2(mask)

    # tuple of axes / negative axes -> let fallback handle
    return NotImplemented


def _build_numpy_function_dispatch():
    import numpy as np
    d = {
        np.concatenate:   _np_concatenate,
        np.vstack:        _np_vstack,
        np.hstack:        _np_hstack,
        np.stack:         _np_stack,
        np.column_stack:  _np_column_stack,
        np.where:         _np_where,
        np.reshape:       _np_reshape,
        np.transpose:     _np_transpose,
        np.dot:           _np_dot,
        np.matmul:        _np_dot,
        np.inner:         _np_inner,
        np.einsum:        _np_einsum,
        np.outer:         _np_outer,
        np.kron:          lambda a, b: kron(a, b),
        np.cross:         _np_cross,
        np.hsplit:        _np_hsplit,
        np.vsplit:        _np_vsplit,
        np.diag:          lambda v, k=0: diag(v) if k == 0 else NotImplemented,
        np.trace:         _np_trace,
        np.repeat:        _np_repeat,
        np.tile:          _np_tile,
        np.zeros_like:    _np_zeros_like,
        np.ones_like:     _np_ones_like,
        np.full_like:     _np_full_like,
        np.cumsum:        _np_cumsum,
        np.sum:           _np_sum,
        np.max:           _np_max,
        np.min:           _np_min,
        np.amax:          _np_max,
        np.amin:          _np_min,
        np.all:           _np_all,
        np.any:           _np_any,
        np.linspace:      _np_linspace,
        np.linalg.norm:   _np_linalg_norm,
        np.linalg.det:    lambda a: det(a),
        np.linalg.inv:    lambda a: inv(a),
        np.linalg.solve:  lambda a, b: solve(a, b),
        np.linalg.cholesky: lambda a: chol(a).T,
        np.clip:          _np_clip,
        np.diff:          _np_diff,
        np.gradient:      _np_gradient,
        np.roll:          _np_roll,
        np.atleast_1d:    _np_atleast_1d,
        np.atleast_2d:    _np_atleast_2d,
        # --- batch 2: synthesized handlers ---
        # --- batch 3: round-3 handlers ---
        # --- batch 4: round-4 handlers ---
        np.pad: _np_pad,
        np.resize: _np_resize,
        np.split: _np_split,
        np.emath.arccos: _np_emath_arccos,
        np.emath.arcsin: _np_emath_arcsin,
        np.emath.arctanh: _np_emath_arctanh,
        np.emath.log2: _np_emath_log2,
        np.emath.logn: _np_emath_logn,
        np.emath.power: _np_emath_power,
        np.choose: _np_choose,
        np.fill_diagonal: _np_fill_diagonal,
        np.select: _np_select,
        np.take: _np_take,
        np.tril_indices_from: _np_tril_indices_from,
        np.triu_indices_from: _np_triu_indices_from,
        np.array_equiv: _np_array_equiv,
        np.iscomplex: _np_iscomplex,
        np.isreal: _np_isreal,
        np.count_nonzero: _np_count_nonzero,
        np.append: _np_append,
        np.array_split: _np_array_split,
        np.asfortranarray: _np_asfortranarray,
        np.block: _np_block,
        np.broadcast_arrays: _np_broadcast_arrays,
        np.broadcast_to: _np_broadcast_to,
        np.delete: _np_delete,
        np.flip: _np_flip,
        np.fliplr: _np_fliplr,
        np.flipud: _np_flipud,
        np.insert: _np_insert,
        np.moveaxis: _np_moveaxis,
        np.ndim: _np_ndim,
        np.ravel: _np_ravel,
        np.rollaxis: _np_rollaxis,
        np.rot90: _np_rot90,
        np.shape: _np_shape,
        np.size: _np_size,
        np.squeeze: _np_squeeze,
        np.swapaxes: _np_swapaxes,
        np.diagonal: _np_diagonal,
        np.linalg.cond: _np_linalg_cond,
        np.linalg.lstsq: _np_linalg_lstsq,
        np.linalg.matrix_power: _np_linalg_matrix_power,
        np.linalg.multi_dot: _np_linalg_multi_dot,
        np.linalg.pinv: _np_linalg_pinv,
        np.linalg.qr: _np_linalg_qr,
        np.linalg.slogdet: _np_linalg_slogdet,
        np.tensordot: _np_tensordot,
        np.vdot: _np_vdot,
        np.allclose: _np_allclose,
        np.array_equal: _np_array_equal,
        np.iscomplexobj: _np_iscomplexobj,
        np.isfortran: _np_isfortran,
        np.isneginf: _np_isneginf,
        np.isposinf: _np_isposinf,
        np.convolve: _np_convolve,
        np.cumprod: _np_cumprod,
        np.interp: _np_interp,
        np.nan_to_num: _np_nan_to_num,
        np.nancumsum: _np_nancumsum,
        np.nanmax: _np_nanmax,
        np.nanmin: _np_nanmin,
        np.real_if_close: _np_real_if_close,
        np.average: _np_average,
        np.corrcoef: _np_corrcoef,
        np.cov: _np_cov,
        np.digitize: _np_digitize,
        np.mean: _np_mean,
        np.ptp: _np_ptp,
        np.std: _np_std,
        np.var: _np_var,
        np.arange: _np_arange,
        np.asanyarray: _np_asanyarray,
        np.asarray: _np_asarray,
        np.ascontiguousarray: _np_ascontiguousarray,
        np.bmat: _np_bmat,
        np.copy: _np_copy,
        np.diagflat: _np_diagflat,
        np.empty_like: _np_empty_like,
        np.full: _np_full,
        np.geomspace: _np_geomspace,
        np.logspace: _np_logspace,
        np.meshgrid: _np_meshgrid,
        np.tril: _np_tril,
        np.triu: _np_triu,
        np.vander: _np_vander,
        np.isclose: _np_isclose,
        np.isrealobj: _np_isrealobj,
        np.isfortran: _np_isfortran,
        np.isscalar: _np_isscalar,
        np.around: _np_around,
        np.ediff1d: _np_ediff1d,
        np.fix: _np_fix,
        np.nanprod: _np_nanprod,
        np.nansum: _np_nansum,
        np.prod: _np_prod,
        np.round: _np_round,
        np.sinc: _np_sinc,
        np.isin: _np_isin,
        # casadi has no complex numbers; np.real/imag/conj are not
        # bridged.  Calls on DM fall through to numpy via .full();
        # symbolic types raise a clear error.
    }
    # numpy.diag dispatches via __array_function__ and ours respects k=0.
    d[np.diag] = _np_diag
    # np.trapz was renamed to np.trapezoid in numpy 2.x and removed from
    # numpy 2.4+.  Reference both names via hasattr so the dict-build
    # doesn't AttributeError on either version: if it did, the
    # `except Exception` in `_numpy_array_function_dispatch` would
    # zero out the dispatch dict and silently break every NEP-18 call.
    for _name in ("trapz", "trapezoid"):
        _f = getattr(np, _name, None)
        if _f is not None:
            d[_f] = _np_trapz
    # np.prod / np.cumprod / np.argmax / np.argmin / np.sort / np.unique / np.eig
    # have no casadi equivalent: numpy fallback (or NotImplemented for symbolic).
    return d


def _numpy_array_function_dispatch(self, func, types, args, kwargs):
    global _NUMPY_FUNCTION_DISPATCH
    if _NUMPY_FUNCTION_DISPATCH is None:
        try:
            _NUMPY_FUNCTION_DISPATCH = _build_numpy_function_dispatch()
        except Exception as e:
            # Building the dispatch dict failed; surface a clear warning
            # rather than silently zeroing the dict (which would make
            # every NEP-18 call on a symbolic type raise the cryptic
            # "no implementation found" TypeError).  Then continue with
            # an empty dict so DM ops still fall back via .full().
            import warnings
            warnings.warn(
                "casadi: failed to build numpy NEP-18 dispatch table "
                "(%s: %s); symbolic numpy ops will be unavailable."
                % (type(e).__name__, e), RuntimeWarning, stacklevel=2)
            _NUMPY_FUNCTION_DISPATCH = {}

    # Gateway (issue #2959): the casadi-aware numpy support dispatches the
    # call through the numpy-semantics array (self wrapped), so the result
    # follows numpy's shape/axis contract and is a NumpyArray.  `self` is
    # the casadi value numpy invoked the protocol on, so it is always a
    # valid proxy even when the casadi operands are nested inside a list
    # argument (np.vstack([M, M]), np.concatenate([...]), ...).
    _mode = GlobalOptions.getNumpyMode()
    if _mode == 1:
        return ArrayInterface._wrap(self, 2).__array_function__(func, types, args, kwargs)

    # Default (legacy casadi 3.7.2) mode: numeric inputs densify to a numpy
    # result; symbolic inputs go through the casadi NEP-18 handler and return
    # a casadi value.  Mode 0 also emits a FutureWarning; mode -1 stays silent.
    if _mode == 0:
        _warn_preserve_type()
    if any(t in (SX, MX) for t in types):
        handler = _NUMPY_FUNCTION_DISPATCH.get(func)
        if handler is not None:
            try:
                result = handler(*args, **kwargs)
            except (TypeError, AttributeError):
                result = NotImplemented
            if result is not NotImplemented:
                return result
        return NotImplemented
    new_args = tuple(a.full() if isinstance(a, DM)
                     else [x.full() if isinstance(x, DM) else x for x in a]
                          if isinstance(a, (list, tuple))
                          and any(isinstance(x, DM) for x in a)
                     else a
                     for a in args)
    return func(*new_args, **kwargs)
