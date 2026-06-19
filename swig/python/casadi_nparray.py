"""Experimental numpy-semantics array wrapper around casadi values.

This module is spliced into the generated casadi module via
`%pythoncode "casadi_nparray.py"` in `swig/casadi.i`, right after the
numpy bridge, so unqualified names (DM/SX/MX, repmat, mtimes, the
`_NUMPY_*_DISPATCH` tables, `_casadi`, ...) resolve against casadi's
own scope.

Motivation (issue #2959): casadi values are *matrices* -- always 2-D,
MATLAB-style single-index semantics.  numpy users expect *arrays* with
a variable number of dimensions, row-major single-index access and row
iteration.  Rather than bend the casadi matrix types (which MATLAB also
depends on), `casadi.ArrayInterface` offers an opt-in, python-only
"array view"; the class doubles as the constructor:

    import casadi as ca
    x = ca.ArrayInterface([[1, 2, 3], [4, 5, 6]])   # logical shape (2, 3)
    x[0]          # -> array of shape (3,)   (numpy row, not a column)
    for row in x: ...                         # iterates rows
    np.sin(x)                                 # -> ArrayInterface
    x @ ca.ArrayInterface([1, 1, 1])          # numpy matmul semantics

Logical ndim 0/1/2 store the data in `_v` as a 2-D casadi value
(sparsity-preserving):
    ndim 0 (scalar)     <-> casadi shape (1, 1)
    ndim 1 (vector n)   <-> casadi shape (1, n)   [a ROW]
    ndim 2 (matrix m,n) <-> casadi shape (m, n)
Storing 1-D as a row makes numpy's right-aligned broadcasting fall out
of casadi's column-broadcast: we repmat operands up to the broadcast
2-D shape ourselves (casadi broadcasts columns but not rows).

DENSE ndim >= 3 is also supported (reshape / transpose(axes) /
broadcasting elementwise ops / axis reductions / stack / concatenate):
`_v` becomes a flat (numel, 1) row-major column and `_shapeN` holds the
logical shape.  An N-D array cannot be handed back to a casadi API
(casadi is 2-D).  Deliberately experimental; API may change.
"""

def _np():
    import numpy
    return numpy


def _is_casadi(x):
    return isinstance(x, (DM, SX, MX))


def _bc_dim(a, b):
    # numpy broadcast of two axis lengths.
    if a == b:
        return a
    if a == 1:
        return b
    if b == 1:
        return a
    raise ValueError("operands could not be broadcast together: %d vs %d" % (a, b))


class ArrayInterface(object):
    """numpy-semantics view over a casadi DM/SX/MX.

    Logical ndim 0/1/2 are stored in `_v` as a 2-D casadi value (the fast,
    sparsity-preserving path).  ndim >= 3 (DENSE only) is stored as a flat
    `_v` column (numel, 1) in row-major / C order plus a `_shapeN` tuple; it
    supports reshape, transpose, broadcasting elementwise ops, axis
    reductions and stack/concatenate, but cannot be handed back to a casadi
    API (casadi is 2-D)."""

    # Win dispatch in mixed binary ops: above plain numpy arrays (priority 0)
    # AND above casadi DM/SX/MX (999/1001/1002), so that `MX + wrapper` defers
    # to our __radd__ (casadi's operators yield to a higher __array_priority__,
    # see _array_priority_yield in casadi.i) and the wrapper wins both ways.
    __array_priority__ = 1005.0

    _shapeN = None   # set (a tuple) only for ndim >= 3

    def __new__(cls, value=None, ndim=None):
        # Dispatch to the backing-typed subclass (ArrayInterfaceDM / SX / MX)
        # so each instance exposes exactly one of __DM__/__SX__/__MX__ -- the
        # one matching its casadi value -- which makes it convert (and
        # type-check) precisely like that casadi type.  `value is None` is the
        # blank-instance path used by _wrap/_wrap_nd, which pick the subclass
        # themselves.
        if value is None or cls is not ArrayInterface:
            return object.__new__(cls)
        if isinstance(value, ArrayInterface):
            v = value._v
        elif _is_casadi(value):
            v = value
        else:
            v = None                      # numpy/list/scalar -> DM-backed
        return object.__new__(_class_for(v))

    def __init__(self, value, ndim=None):
        # numpy constructs by value: copy the casadi backing so that a later
        # in-place __setitem__ never mutates the caller's DM/SX/MX (the
        # internal, copy-free path is _wrap, used for freshly computed values).
        if isinstance(value, ArrayInterface):
            self._v = type(value._v)(value._v)
            self._ndim = value._ndim if ndim is None else ndim
            self._shapeN = value._shapeN
            return
        if _is_casadi(value):
            v = type(value)(value)
            nd = 2 if ndim is None else ndim
        else:
            # python scalar / list / numpy array.  Parse dimensionality
            # with numpy, then move the numbers into a casadi DM.
            np = _np()
            a = np.asarray(value)
            nd = a.ndim if ndim is None else ndim
            if a.ndim == 0:
                v = DM(float(a))
            elif a.ndim == 1:
                v = DM(a.reshape(1, -1))      # store 1-D as a row
            elif a.ndim == 2:
                v = DM(a)
            else:                              # dense N-D (row-major flat)
                self._v = DM(a.reshape(-1, 1, order="C"))
                self._ndim = a.ndim
                self._shapeN = tuple(int(d) for d in a.shape)
                return
        self._v = v
        self._ndim = self._canon_ndim(v, nd)

    # -- construction helpers -------------------------------------------------

    @classmethod
    def _wrap(cls, v, ndim):
        # Build directly from a casadi value already in canonical form for
        # the given ndim (no orientation fix-up beyond _canon); the subclass
        # follows the casadi value's type.
        out = object.__new__(_class_for(v))
        out._v = v
        out._ndim = cls._canon_ndim(v, ndim)
        return out

    # Typed constructors, also the form `repr` produces, e.g.
    # `ArrayInterface.DM([[1, 2], [3, 4]])` or, for a logical shape that
    # differs from the casadi 2-D one, `ArrayInterface.DM([0, 1, 2],
    # shape=(3,))`.  `value` may be anything the casadi type accepts (incl.
    # the casadi value itself).
    @staticmethod
    def _typed(casadi_type, value, shape):
        a = ArrayInterface(value if _is_casadi(value) else casadi_type(value))
        return a if shape is None else _nf_reshape(a, tuple(shape))

    @classmethod
    def DM(cls, value, shape=None):
        return cls._typed(DM, value, shape)

    @classmethod
    def SX(cls, value, shape=None):
        return cls._typed(SX, value, shape)

    @classmethod
    def MX(cls, value, shape=None):
        return cls._typed(MX, value, shape)

    @staticmethod
    def _canon_ndim(v, ndim):
        # Normalise a casadi 2-D value to the canonical storage for `ndim`.
        # (Mutates nothing; callers pass the value back via _wrap/_v.)
        return int(ndim)

    def _canon(self):
        # Ensure self._v matches the canonical orientation for self._ndim.
        if self._ndim == 1 and self._v.size1() != 1:
            self._v = self._v.T          # column -> row
        return self

    # -- basic introspection --------------------------------------------------

    @property
    def ndim(self):
        return len(self._shapeN) if self._shapeN is not None else self._ndim

    @property
    def shape(self):
        if self._shapeN is not None:
            return self._shapeN
        s1, s2 = self._v.size1(), self._v.size2()
        if self._ndim == 0:
            return ()
        if self._ndim == 1:
            return (s2,)
        return (s1, s2)

    # -- N-D (dense, ndim>=3) support -----------------------------------------

    @classmethod
    def _wrap_nd(cls, flat, shape):
        # Build from a flat row-major (C-order) casadi column + logical
        # shape.  Collapses to the fast 0/1/2-D path when rank <= 2.
        shape = tuple(int(d) for d in shape)
        if len(shape) <= 2:
            if len(shape) == 0:
                return cls._wrap(reshape(flat, 1, 1), 0)
            if len(shape) == 1:
                return cls._wrap(reshape(flat, 1, shape[0]), 1)   # (1,n) row
            # (m,n) row-major: flat is C-order -> fill columns of (n,m), T
            return cls._wrap(reshape(flat, shape[1], shape[0]).T, 2)
        v = flat if flat.size2() == 1 else _casadi.vec(flat)
        out = object.__new__(_class_for(v))
        out._v = v
        out._ndim = len(shape)
        out._shapeN = shape
        return out

    def _flat(self):
        # Row-major (C-order) flat column of this array's logical data.
        if self._shapeN is not None:
            return self._v
        nd = self._ndim
        if nd <= 1:
            return _casadi.vec(self._v)            # (1,1)/(1,n) -> column
        return _casadi.vec(self._v.T)              # (m,n) row-major flatten

    def _to_nd(self):
        return self._flat(), self.shape

    @staticmethod
    def _gather(flat, idx):
        # Select elements of a casadi column `flat` at flat indices `idx`,
        # returning a (len(idx), 1) column.  Works for DM/SX/MX.  reshape
        # forces a column even when `flat` is a 1x1 (orientation-ambiguous).
        idx = list(idx)
        if not idx:
            return flat[0:0]
        g = flat[idx]
        return reshape(g, g.numel(), 1)

    @staticmethod
    def _nd_broadcast(arrs):
        # Broadcast a list of arrays to a common shape and return their
        # row-major flat columns (all of length prod(T)) plus T.
        np = _np()
        T = np.broadcast_shapes(*[a.shape for a in arrs])
        n = ArrayInterface._prod(T)
        cols = []
        for a in arrs:
            flat, sh = a._to_nd()
            if sh == T:
                cols.append(flat)
            elif n == 0:
                cols.append(flat)
            else:
                src = np.broadcast_to(
                    np.arange(ArrayInterface._prod(sh)).reshape(sh), T
                ).reshape(-1, order="C").tolist()
                cols.append(ArrayInterface._gather(flat, src))
        return cols, T

    def _nd_reduce(self, axes, mean=False):
        # Sum (or mean) over `axes` of an N-D array; drops those axes.
        np = _np()
        flat, sh = self._to_nd()
        nd = len(sh)
        if axes is None:
            axes = tuple(range(nd))
        elif isinstance(axes, int):
            axes = (axes,)
        axes = tuple(int(x) % nd for x in axes)
        keep = [i for i in range(nd) if i not in axes]
        red_size = ArrayInterface._prod([sh[i] for i in axes])
        keep_shape = tuple(sh[i] for i in keep)
        keep_size = ArrayInterface._prod(keep_shape)
        perm = keep + list(axes)
        src = np.transpose(np.arange(ArrayInterface._prod(sh)).reshape(sh),
                           perm).reshape(-1, order="C").tolist()
        g = ArrayInterface._gather(flat, src)        # C-order: keep..., then red...
        # columns of the (red_size, keep_size) col-major reshape are groups
        s = _casadi.sum1(reshape(g, red_size, keep_size))   # (1, keep_size)
        res = s.T
        if mean and red_size:
            res = res / red_size
        return ArrayInterface._wrap_nd(res, keep_shape)

    @property
    def size(self):
        return self._v.numel()

    @property
    def T(self):
        if self._ndim < 2:
            return self            # numpy: transpose of <2-D is a no-op
        return ArrayInterface._wrap(self._v.T, 2)

    def reshape(self, *shape, **kw):
        if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
            shape = tuple(shape[0])
        return _nf_reshape(self, shape, order=kw.get("order", "C"))

    def transpose(self, *axes):
        if not axes:
            axes = None
        elif len(axes) == 1 and isinstance(axes[0], (tuple, list)):
            axes = tuple(axes[0])
        return _nf_transpose(self, axes)

    def squeeze(self, axis=None):
        # numpy.squeeze: drop length-1 axes (all of them, or the given axis/es).
        shp = self.shape
        if axis is None:
            new = tuple(d for d in shp if d != 1)
        else:
            axs = (axis,) if isinstance(axis, int) else tuple(axis)
            axs = tuple(a % len(shp) for a in axs)
            if any(shp[a] != 1 for a in axs):
                raise ValueError("cannot select an axis to squeeze out which "
                                 "has size not equal to one")
            new = tuple(d for i, d in enumerate(shp) if i not in axs)
        return _nf_reshape(self, new)

    def _natural(self):
        # The underlying casadi value in casadi's vector convention: a
        # logical 1-D (stored as a (1,n) row for numpy broadcasting) maps to
        # a COLUMN (n,1); 0-d and 2-D map to their stored 2-D value.
        if self._shapeN is not None:
            raise TypeError(
                "a %d-D experimental array cannot be used as a casadi value "
                "(casadi is 2-D); reshape to <=2-D first." % self.ndim)
        return self._v.T if self._ndim == 1 else self._v

    def to_casadi(self):
        """The underlying casadi value (1-D -> column (n,1), else 2-D).
        Rarely needed: ca.DM/SX/MX(x) and Function(...) already accept an
        ArrayInterface directly via __DM__/__SX__/__MX__."""
        return self._natural()

    def _hook(self):
        # The casadi value the per-subclass typemap hook (__DM__/__SX__/__MX__)
        # returns: the natural value, or None for an N-D array (which has no
        # 2-D casadi form), so casadi's converter rejects it cleanly.
        return None if self._shapeN is not None else self._natural()

    def to_DM(self):
        return DM(self._v)

    def __len__(self):
        if self._ndim == 0:
            raise TypeError("len() of unsized object")
        return self.shape[0]

    # -- conversion -----------------------------------------------------------

    def __array__(self, dtype=None, copy=None):
        # Densify to numpy (numeric only); enables np.array(x), matplotlib...
        # We always materialise a fresh array, so numpy 2.0's copy=False
        # ("return a view, never copy") cannot be honoured.
        if copy is False:
            raise ValueError("casadi ArrayInterface cannot be converted to a "
                             "numpy view without copying; use copy=True/None")
        if not hasattr(self._v, "full"):
            raise TypeError(
                "cannot convert symbolic %s to a numpy array"
                % type(self._v).__name__)
        np = _np()
        if self._shapeN is not None:          # N-D: flat row-major -> shape
            a = np.array(self._v.full()).reshape(self._shapeN, order="C")
        else:
            a = np.array(self._v.full())      # (m, n)
            if self._ndim == 0:
                a = a.reshape(())
            elif self._ndim == 1:
                a = a.reshape(-1)
        if dtype is not None:
            a = a.astype(dtype)
        return a

    # -- indexing -------------------------------------------------------------

    @staticmethod
    def _norm_int(k, n):
        k = int(k)
        if k < 0:
            k += n
        if k < 0 or k >= n:
            raise IndexError("index %d out of bounds for axis of size %d" % (k, n))
        return k

    def _axis_spec(self, key, n):
        # Translate a numpy axis key into a casadi index + reduces-flag.
        if isinstance(key, slice):
            return key, False
        if isinstance(key, (list, tuple)):
            return [self._norm_int(k, n) for k in key], False
        np = _np()
        if isinstance(key, np.ndarray):
            return [self._norm_int(k, n) for k in key.tolist()], False
        return self._norm_int(key, n), True   # plain int reduces the axis

    def __getitem__(self, idx):
        if self._shapeN is not None:
            raise NotImplementedError(
                "indexing a %d-D experimental array is not supported yet; "
                "reshape to <=2-D first." % self.ndim)
        nd = self._ndim
        keys = idx if isinstance(idx, tuple) else (idx,)
        if nd == 0:
            if keys == ():
                return self
            raise IndexError("too many indices for 0-d array")
        if len(keys) > nd:
            raise IndexError("too many indices: array is %d-d" % nd)

        s1, s2 = self._v.size1(), self._v.size2()
        if nd == 1:
            # stored as (1, n): the single logical axis maps to columns.
            cspec, creduce = self._axis_spec(keys[0], s2)
            sub = self._v[slice(None), cspec]
            res_ndim = 0 if creduce else 1
        else:  # nd == 2
            row_key = keys[0]
            col_key = keys[1] if len(keys) == 2 else slice(None)
            rspec, rreduce = self._axis_spec(row_key, s1)
            cspec, creduce = self._axis_spec(col_key, s2)
            sub = self._v[rspec, cspec]
            res_ndim = (0 if rreduce else 1) + (0 if creduce else 1)
        return ArrayInterface._wrap(sub, res_ndim)._canon()

    def __setitem__(self, idx, value):
        nd = self._ndim
        keys = idx if isinstance(idx, tuple) else (idx,)
        rhs = value._v if isinstance(value, ArrayInterface) else value
        if nd == 0:
            raise IndexError("too many indices for 0-d array")
        s1, s2 = self._v.size1(), self._v.size2()
        if nd == 1:
            cspec, _ = self._axis_spec(keys[0], s2)
            self._v[slice(None), cspec] = rhs
        else:
            row_key = keys[0]
            col_key = keys[1] if len(keys) == 2 else slice(None)
            rspec, _ = self._axis_spec(row_key, s1)
            cspec, _ = self._axis_spec(col_key, s2)
            self._v[rspec, cspec] = rhs

    def __iter__(self):
        if self._ndim == 0:
            raise TypeError("iteration over a 0-d array")
        for i in range(self.shape[0]):
            yield self[i]

    # -- flattening -----------------------------------------------------------

    def _c_order_flat(self):
        # The values in numpy C-order (row-major) as a casadi (numel, 1) column.
        if self._shapeN is not None:
            return self._v                    # N-D storage is already C-order
        return _casadi.vec(self._v.T)         # 2-D: vec(transpose) is row-major

    @property
    def flat(self):
        # numpy.ndarray.flat: a 1-D iterator/accessor over the elements in
        # C-order (supports iteration, len, x.flat[i] and x.flat[i] = v).
        return _FlatIter(self)

    def flatten(self):
        # 1-D copy in C-order (like numpy.ndarray.flatten).
        f = self._c_order_flat()
        return ArrayInterface._wrap(type(f)(f), 1)._canon()

    def ravel(self):
        # 1-D in C-order (like numpy.ravel); a copy here, casadi has no views.
        return self.flatten()

    # -- broadcasting ---------------------------------------------------------

    @staticmethod
    def _as_array(x):
        return x if isinstance(x, ArrayInterface) else ArrayInterface(x)

    @staticmethod
    def _broadcast(arrs):
        # arrs: list of ArrayInterface.  Returns (list_of_casadi_2d, result_ndim)
        # with every value repmat'ed to the common broadcast 2-D shape.
        R = C = 1
        nd = 0
        for a in arrs:
            R = _bc_dim(R, a._v.size1())
            C = _bc_dim(C, a._v.size2())
            nd = max(nd, a._ndim)
        out = []
        for a in arrs:
            v = a._v
            r, c = v.size1(), v.size2()
            if (r, c) != (R, C):
                v = repmat(v, R // r, C // c)
            out.append(v)
        return out, nd

    def _binop(self, other, op):
        a, b = ArrayInterface._as_array(self), ArrayInterface._as_array(other)
        if a.ndim > 2 or b.ndim > 2:
            cols, T = ArrayInterface._nd_broadcast([a, b])
            return ArrayInterface._wrap_nd(op(cols[0], cols[1]), T)
        (va, vb), nd = ArrayInterface._broadcast([a, b])
        return ArrayInterface._wrap(op(va, vb), nd)._canon()

    def _rbinop(self, other, op):
        a, b = ArrayInterface._as_array(other), ArrayInterface._as_array(self)
        if a.ndim > 2 or b.ndim > 2:
            cols, T = ArrayInterface._nd_broadcast([a, b])
            return ArrayInterface._wrap_nd(op(cols[0], cols[1]), T)
        (va, vb), nd = ArrayInterface._broadcast([a, b])
        return ArrayInterface._wrap(op(va, vb), nd)._canon()

    # arithmetic
    def __add__(self, o):      return self._binop(o, lambda x, y: x + y)
    def __radd__(self, o):     return self._rbinop(o, lambda x, y: x + y)
    def __sub__(self, o):      return self._binop(o, lambda x, y: x - y)
    def __rsub__(self, o):     return self._rbinop(o, lambda x, y: x - y)
    def __mul__(self, o):      return self._binop(o, lambda x, y: x * y)
    def __rmul__(self, o):     return self._rbinop(o, lambda x, y: x * y)
    def __truediv__(self, o):  return self._binop(o, lambda x, y: x / y)
    def __rtruediv__(self, o): return self._rbinop(o, lambda x, y: x / y)
    def __floordiv__(self, o): return self._binop(o, lambda x, y: x // y)
    def __rfloordiv__(self, o):return self._rbinop(o, lambda x, y: x // y)
    def __mod__(self, o):      return self._binop(o, lambda x, y: x % y)
    def __rmod__(self, o):     return self._rbinop(o, lambda x, y: x % y)
    def __pow__(self, o):      return self._binop(o, lambda x, y: x ** y)
    def __rpow__(self, o):     return self._rbinop(o, lambda x, y: x ** y)

    # comparisons (elementwise, like numpy)
    def __lt__(self, o):  return self._binop(o, lambda x, y: x < y)
    def __le__(self, o):  return self._binop(o, lambda x, y: x <= y)
    def __gt__(self, o):  return self._binop(o, lambda x, y: x > y)
    def __ge__(self, o):  return self._binop(o, lambda x, y: x >= y)
    def __eq__(self, o):  return self._binop(o, lambda x, y: x == y)
    def __ne__(self, o):  return self._binop(o, lambda x, y: x != y)
    __hash__ = None

    # unary
    def __neg__(self):
        if self._shapeN is not None:
            return ArrayInterface._wrap_nd(-self._v, self._shapeN)
        return ArrayInterface._wrap(-self._v, self._ndim)
    def __pos__(self):  return self
    def __abs__(self):
        if self._shapeN is not None:
            return ArrayInterface._wrap_nd(fabs(self._v), self._shapeN)
        return ArrayInterface._wrap(fabs(self._v), self._ndim)

    # -- matmul (numpy semantics) --------------------------------------------

    def __matmul__(self, other):
        return self._matmul(ArrayInterface._as_array(other))

    def __rmatmul__(self, other):
        return ArrayInterface._as_array(other)._matmul(self)

    def _matmul(self, other):
        a, b = self, other
        # numpy: 1-D on the left is a row, 1-D on the right is a column;
        # the contracted axis is dropped from the result.
        A = a._v                     # ndim1 already stored as (1, n) row
        B = b._v if b._ndim == 2 else b._v.T   # ndim1 -> column for rhs
        prod = mtimes(A, B)
        res_ndim = (1 if a._ndim == 2 else 0) + (1 if b._ndim == 2 else 0)
        return ArrayInterface._wrap(prod, res_ndim)._canon()

    # -- reductions (numpy axis contract: the reduced axis is dropped) -------

    def _axis_norm(self, axis):
        # Normalise a numpy axis for this array; None / 0-d -> None ("all"),
        # else an int axis in [0, ndim).
        if axis is None or self._ndim == 0:
            return None
        ax = int(axis)
        if ax < 0:
            ax += self._ndim
        if ax < 0 or ax >= self._ndim:
            raise ValueError("axis %r out of bounds for %d-d array"
                             % (axis, self._ndim))
        return ax

    def _count(self, axis):
        ax = self._axis_norm(axis)
        return self._v.numel() if ax is None else self.shape[ax]

    def _atleast2d_row(self):
        # numpy.vstack promotes a 1-D (n,) input to a (1, n) row first.
        # ndim-1 is stored as (1, n) already, ndim-0 as (1, 1); relabel.
        return self if self._ndim == 2 else ArrayInterface._wrap(self._v, 2)

    def _reduce_linear(self, axis, red1, red2):
        # Reductions casadi does natively (red1 over rows -> (1,n), red2 over
        # columns -> (m,1)).  Used for sum.
        ax = self._axis_norm(axis)
        v = self._v
        if ax is None:
            return ArrayInterface._wrap(red2(red1(v)), 0)
        if self._ndim == 1:                       # stored (1,n): reduce all
            return ArrayInterface._wrap(red2(v), 0)
        if ax == 0:
            return ArrayInterface._wrap(red1(v), 1)._canon()    # (1,n) -> (n,)
        return ArrayInterface._wrap(red2(v), 1)._canon()        # (m,1) -> (m,)

    def _reduce_fold(self, axis, scalar_reduce):
        # Reductions with no native casadi axis op (max/min/prod/all/any):
        # split into columns/rows and reduce each to a scalar.  Works for
        # DM/SX/MX alike.
        ax = self._axis_norm(axis)
        v = self._v
        if ax is None:
            return ArrayInterface._wrap(scalar_reduce(vec(v)), 0)
        if self._ndim == 1:
            return ArrayInterface._wrap(scalar_reduce(v.T), 0)
        if ax == 0:                               # reduce rows -> (1,n)
            red = hcat([scalar_reduce(c) for c in horzsplit(v)])
            return ArrayInterface._wrap(red, 1)._canon()
        red = vcat([scalar_reduce(r.T) for r in vertsplit(v)])  # cols -> (m,1)
        return ArrayInterface._wrap(red, 1)._canon()

    def sum(self, axis=None):
        if self._shapeN is not None or isinstance(axis, tuple):
            return self._nd_reduce(axis, mean=False)
        return self._reduce_linear(axis, _casadi.sum1, _casadi.sum2)

    def mean(self, axis=None):
        if self._shapeN is not None or isinstance(axis, tuple):
            return self._nd_reduce(axis, mean=True)
        return self.sum(axis) / self._count(axis)

    def prod(self, axis=None):
        import functools as _ft
        return self._reduce_fold(
            axis, lambda c: _ft.reduce(lambda a, b: a * b, vertsplit(c)))

    def max(self, axis=None):
        return self._reduce_fold(axis, mmax)

    def min(self, axis=None):
        return self._reduce_fold(axis, mmin)

    def ptp(self, axis=None):
        return self.max(axis) - self.min(axis)

    def all(self, axis=None):
        return self._reduce_fold(axis, lambda c: mmin(c != 0))

    def any(self, axis=None):
        return self._reduce_fold(axis, lambda c: mmax(c != 0))

    def var(self, axis=None, ddof=0):
        # keepdims-internally so the centred residual broadcasts correctly
        # (numpy's own np.var(a, axis=1) relies on the same).
        ax = self._axis_norm(axis)
        v = self._v
        if ax is None:
            n = v.numel()
            d = v - _casadi.sum1(_casadi.sum2(v)) / n
            return ArrayInterface._wrap(_casadi.sum1(_casadi.sum2(d * d)) / (n - ddof), 0)
        if self._ndim == 1:
            n = v.size2()
            d = v - _casadi.sum2(v) / n
            return ArrayInterface._wrap(_casadi.sum2(d * d) / (n - ddof), 0)
        if ax == 0:
            m = v.size1()
            d = v - repmat(_casadi.sum1(v) / m, m, 1)
            return ArrayInterface._wrap(_casadi.sum1(d * d) / (m - ddof), 1)._canon()
        n = v.size2()
        d = v - repmat(_casadi.sum2(v) / n, 1, n)
        return ArrayInterface._wrap(_casadi.sum2(d * d) / (n - ddof), 1)._canon()

    def std(self, axis=None, ddof=0):
        var = self.var(axis, ddof)
        return ArrayInterface._wrap(sqrt(var._v), var._ndim)

    # ufunc-method -> reduction handler
    _REDUCE_METHODS = None

    @classmethod
    def _reduce_method(cls, name):
        if cls._REDUCE_METHODS is None:
            cls._REDUCE_METHODS = {
                "add": cls.sum, "sum": cls.sum,
                "multiply": cls.prod, "prod": cls.prod,
                "maximum": cls.max, "fmax": cls.max,
                "minimum": cls.min, "fmin": cls.min,
                "logical_and": cls.all, "logical_or": cls.any,
            }
        return cls._REDUCE_METHODS.get(name)

    # -- numpy dispatch -------------------------------------------------------

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # Elementwise ufuncs: broadcast all operands, apply the casadi op,
        # re-wrap with the broadcast ndim.  reduce: axis-aware reduction.
        import numpy as np
        # Data-dependent kwargs we can't bridge (out / where / initial):
        # densify numeric inputs and let numpy compute (returns a plain
        # ndarray); a symbolic input has no such fallback.
        if kwargs.get("out", None) is not None \
                or kwargs.get("where", True) is not True \
                or kwargs.get("initial", np._NoValue) is not np._NoValue:
            if any(ArrayInterface._arg_has_symbolic(x) for x in inputs):
                return NotImplemented
            dens = tuple(ArrayInterface._densify_arg(x) for x in inputs)
            return getattr(ufunc, method)(*dens, **kwargs)
        name = ufunc.__name__
        if method == "__call__":
            op = _NUMPY_UFUNC_DISPATCH.get(name)
            if op is None:
                return NotImplemented
            arrs = [ArrayInterface._as_array(x) for x in inputs]
            if any(a.ndim > 2 for a in arrs):     # dense N-D elementwise
                cols, T = ArrayInterface._nd_broadcast(arrs)
                try:
                    res = op(*cols)
                except (TypeError, ValueError, NotImplementedError):
                    return NotImplemented
                if isinstance(res, tuple):
                    return tuple(ArrayInterface._wrap_nd(r, T) for r in res)
                return ArrayInterface._wrap_nd(res, T)
            vals, nd = ArrayInterface._broadcast(arrs)
            try:
                res = op(*vals)
            except (TypeError, ValueError, NotImplementedError):
                return NotImplemented
            if isinstance(res, tuple):       # multi-output ufunc (modf, ...)
                return tuple(ArrayInterface._wrap(r, nd)._canon() for r in res)
            return ArrayInterface._wrap(res, nd)._canon()
        if method == "reduce":
            a = ArrayInterface._as_array(inputs[0])
            axis = kwargs.get("axis", 0)
            keepdims = bool(kwargs.get("keepdims", False))
            red = ArrayInterface._reduce_method(name)
            if red is not None:
                return a._keepdims(red(a, axis), axis, keepdims)
            # ufuncs with a casadi reduce in the bridge table but no native
            # ndim method (logaddexp -> logsumexp, ...).
            entry = _NUMPY_UFUNC_REDUCE.get(name)
            if entry is not None and entry[0] is not None:
                try:
                    res = entry[0](a._v, 0 if axis is None else int(axis))
                except (TypeError, ValueError, NotImplementedError):
                    return NotImplemented
                nd = 0 if axis is None else max(0, a._ndim - 1)
                return a._keepdims(ArrayInterface._wrap(res, nd)._canon(),
                                   axis, keepdims)
            return NotImplemented
        if method == "accumulate":
            a = ArrayInterface._as_array(inputs[0])
            entry = _NUMPY_UFUNC_REDUCE.get(name)
            if a._shapeN is not None or entry is None or entry[1] is None:
                return NotImplemented
            try:
                return a._accumulate(kwargs.get("axis", 0), entry[1])
            except (TypeError, ValueError, NotImplementedError):
                return NotImplemented
        return NotImplemented

    def _accumulate(self, axis, accum):
        # numpy ufunc.accumulate: cumulative op along one axis, shape kept.
        # `accum(value, casadi_axis)` is the bridge handler (e.g. cumsum).
        ax = self._axis_norm(0 if axis is None else axis)
        if self._ndim == 1:                  # stored (1,n): the lone axis is 1
            return ArrayInterface._wrap(accum(self._v, 1), 1)._canon()
        return ArrayInterface._wrap(accum(self._v, ax), self._ndim)

    def _keepdims(self, result, axis, keepdims):
        # numpy keepdims=True: re-insert the reduced axis as size 1 so the
        # result keeps the input rank.  `result` is the axis-dropped value.
        if not keepdims:
            return result
        nd = self._ndim
        v = result._v
        if nd <= 1:                      # 1-D/scalar -> keep as (1,) / (1,1)
            return ArrayInterface._wrap(v, nd)
        ax = self._axis_norm(axis)       # nd == 2
        if ax == 1:                      # result stored (1,m) row -> (m,1)
            return ArrayInterface._wrap(v.T, 2)
        return ArrayInterface._wrap(v, 2)   # axis None/0 already (1,1)/(1,n)

    @staticmethod
    def _prod(shape):
        n = 1
        for d in shape:
            n *= int(d)
        return n

    def _wrap_generic(self, v, shape=None):
        # Wrap a casadi result from a generic NEP-18 handler.  `shape`, when
        # given, is numpy's exact output shape (inferred by running the
        # numpy function on shape-only dummies) -- trusted only when its
        # element count matches the casadi result, i.e. the casadi handler
        # didn't take a different (1-D vector) interpretation.  Otherwise
        # collapse a vector/scalar casadi shape to its numpy rank: (1,1) ->
        # 0-d, a row/column -> 1-D, else 2-D.
        if shape is not None and ArrayInterface._prod(shape) == v.numel():
            return ArrayInterface._wrap(v, len(shape))._canon()
        s1, s2 = v.size1(), v.size2()
        nd = 0 if s1 * s2 == 1 else (1 if (s1 == 1 or s2 == 1) else 2)
        return ArrayInterface._wrap(v, nd)._canon()

    @staticmethod
    def _dummy_of_shape(shape):
        # A non-degenerate ndarray of `shape` for output-shape inference:
        # ones, but diagonally dominant for a square 2-D matrix so functions
        # like np.linalg.solve don't fail numerically (which would defeat
        # shape inference) on a singular all-zeros/ones matrix.
        np = _np()
        d = np.ones(shape)
        if len(shape) == 2 and shape[0] == shape[1] and shape[0] > 0:
            d = d + np.eye(shape[0]) * shape[0]
        return d

    @staticmethod
    def _shape_dummy(a):
        # A dummy ndarray of an argument's LOGICAL numpy shape (values
        # irrelevant, so symbolic args work too).
        if isinstance(a, ArrayInterface):
            return ArrayInterface._dummy_of_shape(a.shape)
        if isinstance(a, (DM, SX, MX)):
            return ArrayInterface._dummy_of_shape((a.size1(), a.size2()))
        if isinstance(a, (list, tuple)):
            return type(a)(ArrayInterface._shape_dummy(x) for x in a)
        return a

    @staticmethod
    def _infer_shape(func, args, kwargs):
        # numpy's output shape for these argument shapes: a shape tuple, a
        # list of shape tuples (list-returning funcs like hsplit), or None
        # if it can't be determined cheaply (value-dependent / func rejects
        # zeros).
        np = _np()
        try:
            dummies = [ArrayInterface._shape_dummy(a) for a in args]
            ref = func(*dummies, **kwargs)
        except Exception:
            return None
        if isinstance(ref, (list, tuple)):
            return [np.asarray(e).shape for e in ref]
        return np.asarray(ref).shape

    @staticmethod
    def _as_casadi_arg(a):
        # Present an argument as the casadi value an existing NEP-18 handler
        # expects.  A wrapper's ndim-1 becomes a COLUMN (n,1) -- casadi's
        # vector convention and numpy's treatment of 1-D in linear algebra
        # -- not the (1,n) row we store.  Recurses into list/tuple args
        # (np.vstack([...]), np.concatenate([...])); raw casadi and scalars
        # pass through.
        if isinstance(a, ArrayInterface):
            return a._v.T if a._ndim == 1 else a._v
        if isinstance(a, (list, tuple)):
            return type(a)(ArrayInterface._as_casadi_arg(x) for x in a)
        return a

    @staticmethod
    def _arg_has_symbolic(a):
        if isinstance(a, ArrayInterface):
            return isinstance(a._v, (SX, MX))
        if isinstance(a, (SX, MX)):
            return True
        if isinstance(a, (list, tuple)):
            return any(ArrayInterface._arg_has_symbolic(x) for x in a)
        return False

    @staticmethod
    def _densify_arg(a):
        np = _np()
        if isinstance(a, ArrayInterface):
            return np.array(a._v.full())
        if isinstance(a, DM):
            return np.array(a.full())
        if isinstance(a, (list, tuple)):
            return type(a)(ArrayInterface._densify_arg(x) for x in a)
        return a

    def __array_function__(self, func, types, args, kwargs):
        handler = _NPARRAY_FUNCTION_DISPATCH.get(func)
        if handler is not None:
            return handler(*args, **kwargs)
        # Generic fallback: present args in casadi's natural orientation and
        # route DIRECTLY to the casadi NEP-18 handler table -- never through
        # _numpy_array_function_dispatch, whose issue #2959 densify+warn has
        # no place inside this opt-in numpy world.
        unwrapped = tuple(ArrayInterface._as_casadi_arg(a) for a in args)
        tbl = _NUMPY_FUNCTION_DISPATCH
        if tbl is None:
            try:
                tbl = _build_numpy_function_dispatch()
            except Exception:
                tbl = {}
        cas_handler = tbl.get(func)
        if cas_handler is not None:
            try:
                res = cas_handler(*unwrapped, **kwargs)
            except (TypeError, AttributeError):
                res = NotImplemented
            if res is not NotImplemented:
                inferred = ArrayInterface._infer_shape(func, args, kwargs)
                if _is_casadi(res):
                    shape = inferred if not isinstance(inferred, list) else None
                    return self._wrap_generic(res, shape)
                if isinstance(res, (list, tuple)):    # hsplit/vsplit/...
                    shapes = inferred if isinstance(inferred, list) else None
                    return type(res)(
                        self._wrap_generic(x, shapes[i] if shapes else None)
                        if _is_casadi(x) else x
                        for i, x in enumerate(res))
                return res
        # No casadi handler: densify numeric operands and let numpy run
        # (quietly).  Symbolic operands have no numpy fallback.
        if any(ArrayInterface._arg_has_symbolic(a) for a in args):
            return NotImplemented
        return func(*tuple(ArrayInterface._densify_arg(a) for a in args), **kwargs)

    # -- display --------------------------------------------------------------

    def _backing(self):
        return type(self._v).__name__          # 'DM' / 'SX' / 'MX'

    def __str__(self):
        # numpy-style text for numeric data (correct logical shape); the
        # casadi expression for symbolic data.  casadi prints DM/SX with a
        # leading newline -- strip it so our output stays tidy.
        if hasattr(self._v, "full"):           # DM -> exactly like numpy
            return str(self.__array__())
        if self._shapeN is not None:           # symbolic N-D: flat + shape
            return "%s  (flat; logical shape %s)" % (
                str(self._v).strip("\n"), self.shape)
        return str(self._natural()).strip("\n")  # symbolic 0/1/2-D

    def __repr__(self):
        # Round-trippable, e.g. `ArrayInterface.DM([[1, 2], [3, 4]])`; a
        # logical shape that differs from the casadi 2-D one is carried as a
        # `shape=` kwarg (1-D, N-D).  eval(repr(x)) rebuilds x where the
        # contents are literals (numeric, or symbols already in scope).
        backing = self._backing()                       # 'DM'/'SX'/'MX'
        inner = repr(self._v)
        content = (inner[len(backing) + 1:-1].strip()
                   if inner.startswith(backing + "(") else inner)
        shape_kw = ""
        if self.shape != (self._v.size1(), self._v.size2()):
            shape_kw = ", shape=%r" % (self.shape,)
        return "ArrayInterface.%s(%s%s)" % (backing, content, shape_kw)

    def __bool__(self):
        if self._v.numel() != 1:
            raise ValueError(
                "the truth value of an array with more than one element is "
                "ambiguous; use a.any() or a.all()")
        return bool(self._v)

    def _scalar(self):
        if self._v.numel() != 1:
            raise TypeError(
                "only length-1 arrays can be converted to a python scalar")
        return self._v

    def __float__(self):
        return float(self._scalar())

    def __int__(self):
        return int(self._scalar())


# ---- backing-typed subclasses -------------------------------------------
#
# Each carries exactly ONE conversion hook, the one matching its casadi
# value, so it converts (and, via the typed stub, type-checks) precisely
# like that casadi type:
#   - a DM-backed array exposes only __DM__; casadi's to_ptr<SX>/<MX> reach
#     it through their DM path (DM->SX and DM->MX are valid), so a numeric
#     array is accepted everywhere a DM/SX/MX is expected;
#   - an SX-backed array exposes only __SX__ (accepted only where SX fits);
#   - an MX-backed array exposes only __MX__ (accepted only where MX fits).

class ArrayInterfaceDM(ArrayInterface):
    def __DM__(self):
        return self._hook()


class ArrayInterfaceSX(ArrayInterface):
    def __SX__(self):
        return self._hook()


class ArrayInterfaceMX(ArrayInterface):
    def __MX__(self):
        return self._hook()


def _class_for(v):
    if isinstance(v, SX):
        return ArrayInterfaceSX
    if isinstance(v, MX):
        return ArrayInterfaceMX
    return ArrayInterfaceDM            # DM, or about-to-be-DM numeric


class _FlatIter(object):
    """numpy.ndarray.flat: a 1-D iterator/accessor over an ArrayInterface in
    C-order (row-major).  Supports iteration, len(), x.flat[i] / x.flat[s] and
    in-place x.flat[i] = v (which writes back into the owning array)."""

    def __init__(self, owner):
        self._owner = owner
        self._flat = owner._c_order_flat()    # casadi (numel, 1), C-order

    def __len__(self):
        return self._flat.size1()

    def __iter__(self):
        f = self._flat
        for i in range(f.size1()):
            yield ArrayInterface._wrap(f[i, 0], 0)

    def __getitem__(self, k):
        f = self._flat
        n = f.size1()
        if isinstance(k, slice):
            idx = list(range(*k.indices(n)))
            return ArrayInterface._wrap(f[idx, 0], 1)._canon()
        return ArrayInterface._wrap(f[ArrayInterface._norm_int(k, n), 0], 0)

    def _owner_set(self, i, val):
        # Map a C-order flat index back to the owner's casadi cell and write.
        o = self._owner
        v = val._v if isinstance(val, ArrayInterface) else val
        if o._shapeN is not None or o._ndim == 0:
            o._v[i, 0] = v                    # N-D flat storage / 0-d scalar
        elif o._ndim == 1:
            o._v[0, i] = v                    # stored (1, n)
        else:
            o._v[i // o._v.size2(), i % o._v.size2()] = v   # (m, n) row-major

    def __setitem__(self, k, val):
        n = self._flat.size1()
        if isinstance(k, slice):
            idx = list(range(*k.indices(n)))
            vs = val._v if isinstance(val, ArrayInterface) else val
            scalar = not hasattr(vs, "__len__") and getattr(vs, "numel", lambda: 1)() == 1
            for j, i in enumerate(idx):
                self._owner_set(i, vs if scalar else val[j])
        else:
            self._owner_set(ArrayInterface._norm_int(k, n), val)
        self._flat = self._owner._c_order_flat()   # refresh the cached view


# ---- numpy free-function handlers operating in array (ndim-aware) space ----

def _nf_reshape(a, newshape, order='C'):
    a = ArrayInterface._as_array(a)
    if isinstance(newshape, int):
        newshape = (newshape,)
    nd = len(newshape)
    dims = [int(d) for d in newshape]
    total = a.size
    if -1 in dims:
        known = 1
        for d in dims:
            if d != -1:
                known *= d
        dims[dims.index(-1)] = total // known if known else 0
    if order == 'F':                            # column-major (Fortran)
        if a.ndim > 2 or nd > 2:
            raise ValueError("reshape order='F' is only supported for <=2-D")
        flat = _casadi.vec(a._v)
        v = reshape(flat, 1, dims[0]) if nd == 1 else reshape(flat, dims[0], dims[1])
        return ArrayInterface._wrap(v, nd)._canon()
    # 'C' row-major (numpy default): the row-major flat is invariant under a
    # C-order reshape, so just relabel the shape (handles >2-D).
    return ArrayInterface._wrap_nd(a._flat(), tuple(dims))


def _nf_transpose(a, axes=None):
    a = ArrayInterface._as_array(a)
    if a.ndim <= 2 and axes is None:
        return a.T
    np = _np()
    flat, sh = a._to_nd()
    nd = len(sh)
    axes = tuple(range(nd))[::-1] if axes is None else tuple(int(x) % nd for x in axes)
    src = np.transpose(np.arange(ArrayInterface._prod(sh)).reshape(sh),
                       axes).reshape(-1, order="C").tolist()
    g = ArrayInterface._gather(flat, src) if src else flat
    return ArrayInterface._wrap_nd(g, tuple(sh[ax] for ax in axes))


def _nf_concatenate(arrs, axis=0):
    arrs = [ArrayInterface._as_array(x) for x in arrs]
    nd = max(a.ndim for a in arrs)
    if nd > 2:                                  # dense N-D via index gather
        np = _np()
        if axis < 0:
            axis += nd
        big = vcat([a._flat() for a in arrs])
        blocks, off = [], 0
        for a in arrs:
            sh = a.shape
            blocks.append(np.arange(ArrayInterface._prod(sh)).reshape(sh) + off)
            off += ArrayInterface._prod(sh)
        out = np.concatenate(blocks, axis=axis)
        g = ArrayInterface._gather(big, out.reshape(-1, order="C").tolist())
        return ArrayInterface._wrap_nd(g, out.shape)
    if nd <= 1:
        return ArrayInterface._wrap(hcat([a._v for a in arrs]), 1)._canon()
    if axis == 0:
        return ArrayInterface._wrap(vcat([a._v for a in arrs]), 2)._canon()
    if axis in (1, -1):
        return ArrayInterface._wrap(hcat([a._v for a in arrs]), 2)._canon()
    return NotImplemented


def _nf_stack(arrs, axis=0):
    # numpy.stack: join along a NEW axis (may produce an N-D result).
    arrs = [ArrayInterface._as_array(x) for x in arrs]
    np = _np()
    sh = arrs[0].shape
    n = ArrayInterface._prod(sh)
    big = vcat([a._flat() for a in arrs])
    blocks = [np.arange(n).reshape(sh) + i * n for i in range(len(arrs))]
    out = np.stack(blocks, axis=axis)
    g = ArrayInterface._gather(big, out.reshape(-1, order="C").tolist())
    return ArrayInterface._wrap_nd(g, out.shape)


def _nf_dot(a, b, out=None):
    if out is not None:
        return NotImplemented
    return ArrayInterface._as_array(a)._matmul(ArrayInterface._as_array(b))


def _nf_vstack(tup):
    return _nf_concatenate([ArrayInterface._as_array(x)._atleast2d_row() for x in tup],
                           axis=0)


def _nf_hstack(tup):
    arrs = [ArrayInterface._as_array(x) for x in tup]
    nd = max(a._ndim for a in arrs)
    axis = 0 if nd <= 1 else 1
    return _nf_concatenate(arrs, axis=axis)


# Reduction handlers: forward to the ndim-aware methods (axis dropped per
# numpy).  Built as a name->method map so np.sum/np.amax/... share one path.
def _reduction_handler(method_name):
    def handler(a, axis=None, **kw):
        arr = ArrayInterface._as_array(a)
        m = getattr(arr, method_name)
        if method_name in ("var", "std"):
            res = m(axis=axis, ddof=int(kw.get("ddof", 0)))
        else:
            res = m(axis=axis)
        return arr._keepdims(res, axis, bool(kw.get("keepdims", False)))
    return handler


def _build_nparray_function_dispatch():
    np = _np()
    d = {
        np.reshape:     _nf_reshape,
        np.transpose:   _nf_transpose,
        np.concatenate: _nf_concatenate,
        np.stack:       _nf_stack,
        np.vstack:      _nf_vstack,
        np.hstack:      _nf_hstack,
        np.dot:         _nf_dot,
    }
    # axis-aware reductions (numpy name -> ArrayInterface method)
    reducers = {
        "sum": "sum", "mean": "mean", "prod": "prod", "std": "std",
        "var": "var", "amax": "max", "max": "max", "amin": "min",
        "min": "min", "ptp": "ptp", "all": "all", "any": "any",
    }
    for npname, method in reducers.items():
        f = getattr(np, npname, None)
        if f is not None:
            d[f] = _reduction_handler(method)
    return d


_NPARRAY_FUNCTION_DISPATCH = {}
try:
    _NPARRAY_FUNCTION_DISPATCH = _build_nparray_function_dispatch()
except Exception:
    pass


# The class `ArrayInterface` (no underscore) intentionally lands in casadi's
# top-level namespace -> `casadi.ArrayInterface`, which doubles as the
# constructor: `casadi.ArrayInterface(value)` builds a wrapper from a casadi
# value, python scalar, list or numpy array.  We deliberately do NOT expose
# top-level `array` / `asarray` factories -- those names would shadow numpy's
# under `from casadi import *`.
