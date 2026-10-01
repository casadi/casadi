#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
"""Python power-indexing on top of casadi.Struct: nested lists, dicts, callables, prefixes."""
import pickle
from typing import Any
import casadi as ca
import numpy as np

def _is_int(a):
  return isinstance(a, (int, np.integer)) and not isinstance(a, bool)


def _tuple(k):
  return k if isinstance(k, tuple) else (k,)


class Repeater:
  def __init__(self, e):
    self.e = e


def repeated(e):
  """Payload used for every element of a list-valued power index, e.g. s["x",:] = repeated(12)"""
  return Repeater(e)


class NestedDictLiteral:
  """Power index that expands every level into dicts and lists"""


nesteddict = NestedDictLiteral()


def _unpack(payload, i):
  if isinstance(payload, (list, tuple)):
    if i >= len(payload):
      raise Exception("Rhs out of range. Got list index %s but rhs list is only of length %s." % (i, len(payload)))
    return payload[i]
  return payload.e if isinstance(payload, Repeater) else payload


def _apply(p, r):
  try:
    return p(r)
  except NotImplementedError:
    if not isinstance(r, list): raise
    return p(*r)  # casadi's concatenations are variadic


class Delegater:
  """Index into a matrix dimension labeled by a structure, e.g. x["P", index["x"], :]"""
  def __init__(self, arg, name):
    self.arg, self.name = arg, name

  def __str__(self):
    return "%s[%s]" % (self.name, str(self.arg))

  __repr__ = __str__

  def __call__(self, s):
    if not isinstance(s, Structure):
      raise Exception("Cannot use delayed index with a integer shapestruct argument.")
    return s.f[self.arg]


class DelegaterConstructor:
  def __init__(self, name):
    self.name = name

  def __getitem__(self, arg):
    return Delegater(arg, self.name)


index = DelegaterConstructor("index")
indexf = DelegaterConstructor("indexf")


def _sparsity(shape):
  if _is_int(shape): return ca.Sparsity.dense(int(shape), 1)
  if isinstance(shape, (list, tuple)) and len(shape) in (1, 2): return ca.Sparsity.dense(*shape)
  if isinstance(shape, ca.Sparsity): return shape
  raise Exception("The 'shape' argument, if present, must be an integer, a tuple of 1 or 2 integers, or a sparsity pattern. Got %s " % str(shape))


class StructEntry:
  keys = ['repeat', 'shape', 'sym', 'expr', 'struct', 'shapestruct', 'type']
  conflicts = [('shape', ['struct']), ('struct', ['shape', 'shapestruct']), ('shapestruct', ['struct']),
               ('sym', ['shape', 'repeat', 'expr']), ('expr', ['shape', 'repeat', 'sym'])]

  def __init__(self, *args, **kwargs):
    if len(args) != 1:
      raise Exception("Expected the entry name as only positional argument, got %s" % str(args))
    self.name, self.kwargs = args[0], kwargs
    for k in kwargs:
      if k not in self.keys: raise Exception("Unknown keyword argument '%s'. Please use one of %s." % (k, str(self.keys)))
    for k, fk in self.conflicts:
      for f in fk:
        if k in kwargs and f in kwargs:
          raise Exception("You supplied keyword argument '%s', but it cannot be combined with keyword argument '%s'." % (k, f))
    self.repeat = kwargs.get("repeat", [])
    if not isinstance(self.repeat, list): self.repeat = [self.repeat]
    if not all(_is_int(x) for x in self.repeat):
      raise Exception("The 'repeat' argument, if present, must be a list of integers, but got %s" % str(self.repeat))
    self.repeat = [int(x) for x in self.repeat]
    self.struct, self.sym, self.expr = kwargs.get("struct"), kwargs.get("sym"), kwargs.get("expr")
    sp = _sparsity(kwargs.get("shape", 1))
    self.shapestruct = kwargs.get("shapestruct")
    if self.shapestruct is not None:
      ss = self.shapestruct if isinstance(self.shapestruct, tuple) else (self.shapestruct,)
      if not 0 < len(ss) <= 2 or not all(isinstance(e, Structure) or _is_int(e) for e in ss):
        raise Exception("The 'shapestruct' argument, if present, must be a structure or a tuple of structures or numbers")
      self.shapestruct = ss + (1,) * (2 - len(ss))
    if self.sym is not None:
      if isinstance(self.sym, Structure):
        self.struct, self.sym = self.sym, self.sym.cat
      elif not (isinstance(self.sym, ca.SX) and self.sym.is_valid_input()):
        raise Exception("The 'sym' argument must be a purely symbolic SX or a structured symbolic. Got %s instead." % str(self.sym))
      sp = self.sym.sparsity()
    if self.expr is not None:
      e, self.repeat = self.expr, []
      while isinstance(e, list):
        self.repeat.append(len(e))
        e = e[0] if e else None
      if e is None:
        sp = ca.Sparsity(0, 0)
      elif hasattr(e, "sparsity"):
        sp = e.sparsity()
      else:
        raise Exception("The 'expr' argument must be a matrix expression or nested list of matrix expressions. Got %s instead." % str(e))
    self.type = kwargs.get("type")
    if self.type not in (None, 'symm'):
      raise Exception("You supplied a type argument '%s' but it is not recognised. Use one of ['symm']" % str(self.type))
    symm = self.type == 'symm'
    if symm and sp.size1() != sp.size2():
      raise Exception("You supplied a type 'symm', but matrix is not square. Got %s." % sp.dim())
    if self.struct is not None:
      self.content = self.struct._s
    elif self.shapestruct is not None:
      labels = [e._s if isinstance(e, Structure) else ca.Struct.leaf(ca.Sparsity.dense(int(e), 1)) for e in self.shapestruct]
      self.content = ca.Struct.matrix(labels[0], labels[1], symm)
    else:
      self.content = ca.Struct.leaf(sp, symm)

  def __reduce__(self):
    return (_entry, (self.name, self.kwargs))


def _entry(name, kwargs):
  return StructEntry(name, **kwargs)


def entry(*args, **kwargs):
  if len(args) == 1 and not kwargs and isinstance(args[0], StructEntry): return args[0]
  return StructEntry(*args, **kwargs)


class Structure(object):
  """Layout of a flat vector, optionally holding data: power indexing with s[...]"""

  description = "Generic Structured object"

  def __init__(self, arg, order=None):
    self._v, self._target, self._args = None, None, (arg, order)
    if isinstance(arg, Structure):
      self._s, self._entries, self._declared = arg._s, arg._entries, []
      return
    if not isinstance(arg, list):
      raise Exception("Expecting list of entries, with possible tuples for grouping, but got %s" % str(arg))
    entries, groups = [], []
    for e in arg:
      es = [entry(x) for x in (e if isinstance(e, tuple) else (e,))]
      entries += es
      groups.append(tuple(x.name for x in es) if isinstance(e, tuple) else es[0].name)
    names = [e.name for e in entries]
    duplicates = sorted(set(n for n in names if names.count(n) > 1))
    if duplicates: raise Exception("Your list of entries contains duplicates: %s" % str(duplicates))
    if order is not None:
      if any(isinstance(g, tuple) for g in groups):
        raise Exception("You supplied an order by using tuple syntax on entries, but you overwrite it with the 'order' keyword. Use one or the other, not both.")
      groups = order
    self._entries, self._declared, self._s = dict(zip(names, entries)), entries, ca.Struct()
    for g in groups + [n for n in names if n not in sum([list(_tuple(g)) for g in groups], [])]:
      for n in _tuple(g):
        if n not in self._entries: raise Exception("Order '%s' is invalid." % n)
        self._s.add(n, self._entries[n].content, self._entries[n].repeat)
      if isinstance(g, tuple) and len(g) > 1: self._s.interleave(list(g))

  def _bind(self, v):
    self._v, self._mtype = v, type(v.cat())
    return self

  def __DM__(self) -> Any:
    return _cast(self.cat, ca.DM) if self._v else None

  def __SX__(self) -> Any:
    return _cast(self.cat, ca.SX) if self._v else None

  def __MX__(self) -> Any:
    return _cast(self.cat, ca.MX) if self._v else None

  def _entry(self, name):
    if name not in self._entries:
      raise Exception("Unknown keyword '%s'. Available entries: %s" % (name, str(self.keys())))
    return self._entries[name]

  struct = property(lambda self: self)
  size = property(lambda self: self._s.nnz())
  shape = property(lambda self: (self.size, 1))
  cat = property(lambda self: self._v.cat() if self._target is None else self._target)
  prefix = property(lambda self: PrefixConstructor(self))
  i = property(lambda self: _IndexGetter(self, False))
  f = property(lambda self: _IndexGetter(self, True))
  map = property(lambda self: _IndexMap(self))

  def sparsity(self):
    return ca.Sparsity.dense(self.size, 1)

  def keys(self):
    return list(self._s.names())

  def getCanonicalIndex(self, i, extraMode=1):
    """Entry path of flat element i, plus its nonzero index (1) or (column, row) (2) in the matrix"""
    p, n, s = list(self._s.path(i)), 0, self
    while s is not None:
      e = s._entry(p[n])
      n, s = n + 1 + len(e.repeat), e.struct
    can = tuple(p[:n])
    if extraMode == 0: return can
    k = i - int(min(self._s.index(list(can)).nonzeros()))
    if extraMode == 1: return can + (k,)
    sp = ca.Sparsity.triu(e.content.sparsity()) if e.type == "symm" else e.content.sparsity()
    return can + (sp.get_col()[k], sp.row()[k])

  def canonicalIndices(self, extraMode=1):
    return [self.getCanonicalIndex(i, extraMode) for i in range(self.size)]

  def getLabel(self, i, extraMode=1):
    return "[" + ",".join(map(str, self.getCanonicalIndex(i, extraMode))) + "]"

  def labels(self, extraMode=1):
    return [self.getLabel(i, extraMode) for i in range(self.size)]

  def __str__(self, compact=False):
    return str(self._s) if self._v is None else str(self._v)

  __repr__ = __str__

  def __reduce__(self):
    return (type(self), self._args)

  def save(self, filename):
    with open(filename, "wb") as f: pickle.dump(self, f, 2)

  def __getitem__(self, pi):
    return self._walk(self, _tuple(pi), [], None, "get")

  def __setitem__(self, pi, value):
    self._walk(self, _tuple(pi), [], value, "set")
    if self._target is not None: self._target[:, :] = ca.reshape(self._v.cat(), self._target.shape)

  def _leaf(self, path, payload, mode):
    if mode == "index": return self._s.index(path)
    if mode == "get": return self._v.get(path)
    self._v.set(path, self._mtype(payload.e if isinstance(payload, Repeater) else payload))

  def _call(self, p, walk, payload, mode):
    if mode != "set": return _apply(p, walk(mode))
    self._v.set_nz(_apply(p, walk("index")), self._mtype(payload))

  def _walk(self, s, pi, path, payload, mode, e=None):
    """Walk power index pi from record s (None at a matrix of entry e), resolving lists, dicts and callables"""
    try:
      if pi and callable(pi[0]) and not isinstance(pi[0], (Structure, Delegater)):
        return self._call(pi[0], lambda m: self._walk(s, pi[1:], path, payload, m, e), payload, mode)
      if s is None or not pi:
        extra = [(p(e.shapestruct[k]) if isinstance(p, Delegater) else p) for k, p in enumerate(pi)
                 if not isinstance(p, NestedDictLiteral)]
        return self._leaf(path + ca.StructDM._path(tuple(extra)), payload, mode)
      p, rest = pi[0], pi[1:]
      if isinstance(p, str):
        return self._repeat(s._entry(p), rest, path + [p], payload, mode, s._entry(p).repeat)
      if p is Ellipsis or isinstance(p, list):
        return [self._repeat(s._entry(k), rest, path + [k], _unpack(payload, i), mode, s._entry(k).repeat)
                for i, k in enumerate(s.keys() if p is Ellipsis else p)]
      if isinstance(p, (dict, set, NestedDictLiteral)):
        r = pi if isinstance(p, NestedDictLiteral) else rest
        keys = [k for k in payload if not isinstance(p, set) or k in p] if isinstance(payload, dict) else (p if isinstance(p, set) else s.keys())
        return dict((k, self._repeat(s._entry(k), r, path + [k], payload[k] if isinstance(payload, dict) else payload, mode, s._entry(k).repeat)) for k in keys)
      if isinstance(p, slice): raise Exception("slice not allowed here, did you mean '...' ?")
      raise Exception("I don't know what to do with this: %s" % str(p))
    except Exception as err:
      raise Exception("Error occured in struct context with powerIndex %s, at canonicalIndex %s" % (str(pi), str(path))) from err

  def _repeat(self, e, pi, path, payload, mode, dims):
    if not dims: return self._walk(e.struct, pi, path, payload, mode, e)
    if not pi: pi = (slice(None),)
    p, rest, n = pi[0], pi[1:], dims[0]
    if callable(p) and not isinstance(p, Structure):
      return self._call(p, lambda m: self._repeat(e, rest, path, payload, m, dims), payload, mode)
    if isinstance(p, slice): p = list(range(*p.indices(n)))
    if _is_int(p): return self._repeat(e, rest, path + [int(p)], payload, mode, dims[1:])
    if isinstance(p, (list, NestedDictLiteral)):
      r = pi if isinstance(p, NestedDictLiteral) else rest
      return [self._repeat(e, r, path + [int(i)], _unpack(payload, j), mode, dims[1:])
              for j, i in enumerate(range(n) if isinstance(p, NestedDictLiteral) else p)]
    if isinstance(p, (dict, set)): raise Exception("powerIndex entry %s cannot be used in list context." % str(p))
    raise Exception("I don't know what to do with this: %s" % str(p))

  def __call__(self, arg: Any = 0) -> Any:
    a, mtype = _argtype(arg)
    return {ca.DM: DMStruct, ca.SX: SXStruct, ca.MX: MXStruct}[mtype](self, data=a)

  def _view(self, e, arg: Any, rows, cols) -> Any:
    a, mtype = _argtype(arg)
    if a.is_scalar() and a.numel() != rows * cols: a = ca.DM.ones(rows, cols) * a
    if a.shape != (rows, cols): raise Exception("Expecting %s of shape (%d,%d). Got %s" % (mtype.__name__, rows, cols, a.dim()))
    s = Structure([e])
    s._bind(getattr(ca, "Struct" + type(a).__name__)(s._s, ca.vec(a)))
    s._target = a
    return Prefixer(s, ("t",), castmaster=True)

  def repeated(self, arg: Any = 0) -> Any:
    a, _ = _argtype(arg)
    return self._view(entry("t", struct=self, repeat=a.shape[1]), a, self.size, a.shape[1])

  def squared(self, arg: Any = 0) -> Any:
    return self._view(entry("t", shapestruct=(self, self)), arg, self.size, self.size)

  def product(self, otherstruct, arg: Any = 0) -> Any:
    return self._view(entry("t", shapestruct=(self, otherstruct)), arg, self.size, otherstruct.size)

  def squared_repeated(self, arg: Any = 0) -> Any:
    a, _ = _argtype(arg)
    return self._view(entry("t", shapestruct=(self, self), repeat=a.shape[1] // max(self.size, 1)), a, self.size, a.shape[1])



def _cast(x, t):
  return x if isinstance(x, t) else None


def _argtype(arg):
  for t in (ca.DM, ca.MX, ca.SX):
    if isinstance(arg, t): return arg, t
  for t in (ca.DM, ca.MX, ca.SX):
    try:
      return t(arg), t
    except Exception:
      pass
  raise Exception("Call to Structure has weird argument: expecting DM-like or MX-like or SXMatrix-like")


class _IndexGetter:
  def __init__(self, s, flat):
    self.s, self.flat = s, flat

  def __getitem__(self, pi):
    r = self.s._walk(self.s, _tuple(pi), [], None, "index")
    return _flat(r) if self.flat else r


def _flat(r):
  if isinstance(r, list): return sum([_flat(e) for e in r], [])
  return [int(i) for i in r.nonzeros()]


class _IndexMap:
  def __init__(self, s):
    self.s = s

  def __getitem__(self, k):
    return self.s._s.index(list(k))


class Prefixer:
  def __init__(self, struct, prefix, castmaster=False):
    self.struct, self.prefix, self.castmaster = struct, prefix, castmaster

  def __DM__(self) -> Any:
    return _cast(self.cast(), ca.DM)

  def __SX__(self) -> Any:
    return _cast(self.cast(), ca.SX)

  def __MX__(self) -> Any:
    return _cast(self.cast(), ca.MX)

  def __reduce__(self):
    return (Prefixer, (self.struct, self.prefix, self.castmaster))

  def __getattr__(self, name):
    # Delegate to the expression this prefix resolves to, e.g. x.shape
    if name in ("struct", "prefix", "castmaster"): raise AttributeError(name)
    t = self.cast()
    if isinstance(t, (list, dict, tuple)):
      raise AttributeError("Cannot get attribute '%s': prefix %s resolves to a %s, not to a matrix." % (name, str(self.prefix), type(t).__name__))
    return getattr(t, name)

  def cast(self):
    return self.struct.cat if self.castmaster else self()

  def __str__(self):
    return "prefix(" + str(self.prefix) + "," + str(self.struct._s) + ")"

  __repr__ = __str__

  def __call__(self):
    return self.struct[self.prefix]

  def __getitem__(self, pi):
    return self.struct[self.prefix + _tuple(pi)]

  def __setitem__(self, pi, data):
    self.struct[self.prefix + _tuple(pi)] = data


class PrefixConstructor:
  def __init__(self, struct):
    self.struct = struct

  def __str__(self):
    return "prefixConstructor(" + str(self.struct._s) + ")"

  __repr__ = __str__

  def __getitem__(self, prefix):
    return Prefixer(self.struct, _tuple(prefix))


class CasadiStructured(Structure):
  """Layout only, e.g. struct(["x", "y"]); call it for numeric or symbolic values"""


class ssymStruct(Structure):
  description = "symbolic SX"

  def __init__(self, arg, order=None):
    Structure.__init__(self, arg, order)
    if any(e.expr is not None for e in self._declared):
      raise Exception("struct_symSX does not accept entries with an 'expr' argument, because such an element is not purely symbolic.")
    self._bind(ca.StructSX.sym(self._s))
    for e in self._declared:
      if e.sym is not None: self._v.set([e.name], e.sym)

  def __setitem__(self, pi, value):
    raise TypeError("'%s' object does not support item assignment" % type(self).__name__)


class msymStruct(ssymStruct):
  description = "MX.sym"

  def __init__(self, arg, order=None):
    Structure.__init__(self, arg, order)
    if any(e.expr is not None for e in self._declared):
      raise Exception("struct_symMX does not accept entries with an 'expr' argument, because such an element is not purely symbolic.")
    if any(e.sym is not None for e in self._declared):
      raise Exception("struct_symMX does not accept entries with an 'sym' argument.")
    self._bind(ca.StructMX.sym(self._s))


class MatrixStruct(Structure):
  """Mutable values: entries given by expr=, or data wrapped; unset parts are NaN"""
  mtype = ca.DM

  def __init__(self, arg, data=None, order=None):
    Structure.__init__(self, arg, order)
    if data is None and any(e.expr is None for e in self._declared):
      raise Exception("struct_%s does only accept entries with an 'expr' argument." % self.mtype.__name__)
    self._args = (arg, data, order)
    self._bind(getattr(ca, "Struct" + self.mtype.__name__)(self._s, self.mtype.nan(1, 1) if data is None else self.mtype(data)))
    for e in self._declared:
      if e.expr is not None: self[e.name] = e.expr

  @property
  def description(self):
    return "Mutable " + self.mtype.__name__

  def __reduce__(self):
    # Numeric values travel along; symbolic ones would need a pickle context
    return (type(self), (self._args[0], self.cat if self.mtype is ca.DM else None, self._args[2]))


class DMStruct(MatrixStruct):
  mtype = ca.DM


class SXStruct(MatrixStruct):
  mtype = ca.SX


class MXStruct(MatrixStruct):
  mtype = ca.MX


class MXVeccatStruct(MatrixStruct):
  mtype = ca.MX

  def __init__(self, arg, order=None):
    MatrixStruct.__init__(self, arg, None, order)
    self._args = (arg, order)

  def __reduce__(self):
    return (type(self), self._args)


struct = CasadiStructured
struct_symSX = ssymStruct
struct_symMX = msymStruct
struct_SX = SXStruct
struct_MX = MXVeccatStruct
struct_MX_mutable = MXStruct


def struct_load(filename):
  with open(filename, "rb") as f: return pickle.load(f)
