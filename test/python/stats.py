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
import casadi as ca
import numpy as np
import unittest
import ctypes
import hashlib
import json
import itertools
import os
import sys
import struct
import subprocess
import tempfile
import time
from helpers import openmp_flags, casadiTestCase, requires_integrator, requires_nlpsol, requiresPlugin, slow


# Minimal CBOR encoder for the records of runtime/casadi_stats.hpp
def _head(major, v):
  if v < 24: return bytes([major * 32 + v])
  for ai, n in ((24, 1), (25, 2), (26, 4), (27, 8)):
    if v < 256 ** n: return bytes([major * 32 + ai]) + v.to_bytes(n, "big")
  raise ValueError(v)

def _cbor(x):
  if x is None: return b"\xf6"
  if x is True: return b"\xf5"
  if x is False: return b"\xf4"
  if isinstance(x, int): return _head(0, x) if x >= 0 else _head(1, -1 - x)
  if isinstance(x, float): return b"\xfb" + struct.pack(">d", x)
  if isinstance(x, str):
    b = x.encode()
    return _head(3, len(b)) + b
  if isinstance(x, list): return _head(4, len(x)) + b"".join(_cbor(e) for e in x)
  raise TypeError(x)

# Record kinds of runtime/casadi_stats.hpp
BEGIN_CALL, END_CALL, SET, APPEND, BEGIN_SECTION, END_SCOPE, BEGIN_ITERATION, \
  DECLARE_FIELDS, SET_FIELD = range(1, 10)

class _Stream:
  def __init__(self):
    self.buf = b""
    self.bounds = {0}
  def rec(self, *r):
    pos = len(self.buf)
    self.buf += _cbor(list(r))
    self.bounds.add(len(self.buf))
    return pos


def sections(d):
  """The calls below a call with their section: pre, (iterations, k), post or children"""
  r = [("pre", c) for c in d.get("pre", {}).get("children", [])]
  r += [(("iterations", k), c) for k, it in enumerate(d.get("iterations", []))
        for c in it["children"]]
  r += [("post", c) for c in d.get("post", {}).get("children", [])]
  return r + [("children", c) for c in d.get("children", [])]

def calls(d):
  return [c for _, c in sections(d)]

def names(l):
  """Names of calls, or of the calls of a section"""
  return [c["name"] for c in (l["children"] if isinstance(l, dict) else l)]

def walk(d):
  yield d
  for c in calls(d): yield from walk(c)

def find(roots, name):
  return [c for r in roots for c in walk(r) if c["name"] == name]

def shape(d):
  """Structure of a call tree: names, sections and flags, no values"""
  return (d["name"], d["flag"], tuple((k, shape(c)) for k, c in sections(d)))

def is_prefix(cut, full):
  """A tree decoded from a cut stream holds a prefix of the full tree"""
  a, b = sections(cut), sections(full)
  if cut["name"] != full["name"] or len(a) > len(b): return False
  if cut["flag"] is not None and cut["flag"] != full["flag"]: return False
  for k, v in cut["stats"].items():
    if k == "iterations":
      # A value not recorded yet (None, NaN) matches whatever the full stream holds
      same = lambda a, b: a is None or a != a or a == b
      if any(not all(map(same, e, full["stats"][k][c])) for c, e in v.items()): return False
    elif full["stats"].get(k) != v:
      return False
  return all(ka == kb and is_prefix(ca_, cb) for (ka, ca_), (kb, cb) in zip(a, b))

def strip_ids(d):
  """A call tree without ids, which depend on the instance"""
  r = {k: v for k, v in d.items() if k != "id"}
  if "children" in d: r["children"] = [strip_ids(c) for c in d["children"]]
  for k in ["pre", "post"]:
    if k in d: r[k] = strip_ids(d[k])
  if "iterations" in d: r["iterations"] = [strip_ids(it) for it in d["iterations"]]
  return r


class _Rosenbrock:
  def __init__(self):
    x = ca.MX.sym('x', 2)
    p = ca.MX.sym('p')
    nlp = {'x': x, 'p': p, 'f': (1-x[0])**2 + p*(x[1]-x[0]**2)**2, 'g': x[0]+x[1]}
    self.solver = ca.nlpsol('solver', 'ipopt', nlp,
      {'print_time': False, 'ipopt.print_level': 0, 'ipopt.sb': 'yes'})
    self.args = dict(x0=[-1, 1], lbg=-ca.inf, ubg=1)


# unsigned char* (*reserve)(struct casadi_stats_sink* s, casadi_int len, casadi_int* pos)
CReserve = ctypes.CFUNCTYPE(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_longlong,
                            ctypes.POINTER(ctypes.c_longlong))

# void (*commit)(struct casadi_stats_sink* s, casadi_int pos)
CCommit = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_longlong)

class CStats(ctypes.Structure):
  # struct casadi_stats_sink
  _fields_ = [("p", ctypes.POINTER(ctypes.c_ubyte)), ("cap", ctypes.c_longlong),
              ("needed", ctypes.c_longlong), ("reserve", CReserve), ("commit", CCommit),
              ("data", ctypes.c_void_p)]

def get_stat(lib, stream, fname, key, kind="int", cap=64):
  """<prefix>_get_stat_int/real/text on stream, or 'missing'"""
  f = getattr(lib, lib._prefix + "_get_stat_" + kind)
  n = len(stream)
  data = (ctypes.c_ubyte * max(n, 1)).from_buffer_copy(stream + b"\0")
  sink = CStats(ctypes.cast(data, ctypes.POINTER(ctypes.c_ubyte)), n, n, CReserve(), CCommit(), None)
  args = [ctypes.byref(sink), fname.encode(), key.encode()]
  if kind == "text":
    v = ctypes.create_string_buffer(b"#" * max(cap, 1), max(cap, 1))
    if f(*args, v, ctypes.c_longlong(cap)): return "missing"
    return v.value.decode() if cap > 0 else v.raw.decode()
  v = ctypes.c_longlong() if kind == "int" else ctypes.c_double()
  if f(*args, ctypes.byref(v)): return "missing"
  return v.value

ROOT = ctypes.c_longlong(-1)

def parallel_options(par):
  """extra_options and definitions for a map of kind par"""
  if par == "openmp":
    return openmp_flags(), ["CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_OMP"]
  if par == "thread":
    if os.name == "nt": return [], ["CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_WINDOWS"]
    return ["-pthread"], ["CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_POSIX"]
  return [], []

def cdll(libname, prefix):
  lib = ctypes.CDLL(os.path.abspath(libname))
  setattr(lib, "_prefix", prefix)
  return lib


class Statstests(casadiTestCase):

  # --- helpers for generated code ---

  def build(self, F, inputs, opts=None, std="c89", extra_options=None, extralibs=None,
            extralibdirs=None, definitions=None, digits=15, main=False, valgrind=False,
            helgrind=False):
    """check_codegen with stats; returns the library loaded with ctypes"""
    o = {"stats": True, "with_header": True}
    o.update(opts or {})
    ret = self.check_codegen(F, inputs=inputs, opts=o, std=std, extra_options=extra_options,
                             extralibs=extralibs or "", extralibdirs=extralibdirs or [],
                             definitions=definitions, digits=digits, main=main,
                             valgrind=valgrind, helgrind=helgrind)
    libname = ret["libname"]
    return cdll(libname, os.path.splitext(os.path.basename(libname))[0])

  def call(self, lib, F, inputs, cap=1 << 16, with_stats=True, reserve=None, commit=None,
           raw=False, unrolled=False):
    """Call F with a sink; reserve/commit: owner callbacks; raw: hand-made sink without them

    unrolled: through <name>_unrolled_with_stats (codegen option unroll_args)

    Returns flag, outputs, StatsRecorder, needed"""
    name = F.name()
    ci = ctypes.c_longlong
    sz = [ci() for _ in range(4)]
    getattr(lib, name + "_work")(*[ctypes.byref(s) for s in sz])
    sz_arg, sz_res, sz_iw, sz_w = [s.value for s in sz]
    ins = [np.ascontiguousarray(np.broadcast_to(np.array(ca.DM(x).nonzeros(), dtype=float),
           F.nnz_in(i))) for i, x in enumerate(inputs)]
    outs = [np.zeros(F.nnz_out(i)) for i in range(F.n_out())]
    P = ctypes.POINTER(ctypes.c_double)
    arg = (P * max(sz_arg, 1))(*[x.ctypes.data_as(P) for x in ins])
    res = (P * max(sz_res, 1))(*[x.ctypes.data_as(P) for x in outs])
    iw = (ci * max(sz_iw, 1))()
    w = (ctypes.c_double * max(sz_w, 1))()
    getattr(lib, name + "_incref")()
    mem = getattr(lib, name + "_checkout")()
    buf = (ctypes.c_ubyte * max(cap, 1))()
    p = ctypes.cast(buf, ctypes.POINTER(ctypes.c_ubyte))
    if raw:
      st = b = CStats(p, cap, 0, CReserve(), CCommit(), None)
    else:
      make = getattr(lib, getattr(lib, "_prefix") + "_stats_make_sink")
      make.restype = CStats
      st = b = make(p, ctypes.c_longlong(cap))
      if reserve: st.reserve = reserve
      if commit: st.commit = commit
    if with_stats and unrolled:
      flag = getattr(lib, name + "_unrolled_with_stats")(
        *[x.ctypes.data_as(P) for x in ins + outs], arg, res, iw, w, mem, ctypes.byref(st), ROOT)
    elif with_stats:
      flag = getattr(lib, name + "_with_stats")(arg, res, iw, w, mem, ctypes.byref(st),
                                                ROOT)
    else:
      flag = getattr(lib, name)(arg, res, iw, w, mem)
    getattr(lib, name + "_release")(mem)
    getattr(lib, name + "_decref")()
    self.stream = bytes(buf)[:min(cap, b.needed)]
    return flag, outs, ca.StatsRecorder.from_bytes(self.stream), b.needed

  # --- decoding ---

  def test_decode_interleaved(self):
    # Interleaved records of two concurrent map iterations
    s = _Stream()
    root = s.rec(BEGIN_CALL, None, "p_f0:G", None)
    m = s.rec(BEGIN_CALL, root, "p_f1:map", None)
    ia = s.rec(BEGIN_ITERATION, m, 1)
    ib = s.rec(BEGIN_ITERATION, m, 0)
    fa = s.rec(BEGIN_CALL, ia, "p_f2:solver", 0)
    fb = s.rec(BEGIN_CALL, ib, "p_f2:solver", 1)
    s.rec(DECLARE_FIELDS, fb, ["obj", "mu"])
    s.rec(DECLARE_FIELDS, fa, ["obj", "mu"])
    a0 = s.rec(BEGIN_ITERATION, fa, 0)
    s.rec(SET_FIELD, a0, 0, 1.5)
    b0 = s.rec(BEGIN_ITERATION, fb, 0)
    s.rec(SET_FIELD, a0, 1, 0.1)
    # Fields in any order
    s.rec(SET_FIELD, b0, 1, 0.2)
    s.rec(SET_FIELD, b0, 0, 2.5)
    s.rec(END_SCOPE, a0)
    a1 = s.rec(BEGIN_ITERATION, fa, 1)
    s.rec(SET_FIELD, a1, 0, 1.0)
    s.rec(SET_FIELD, a1, 1, 0.01)
    s.rec(APPEND, fa, "warnings", "restoration")
    s.rec(SET, fb, "return_status", "Solve_Succeeded")
    s.rec(END_SCOPE, a1)
    s.rec(END_SCOPE, b0)
    s.rec(SET, fa, "iter_count", 2)
    s.rec(SET, fa, "success", True)
    s.rec(END_CALL, fb, 0)
    s.rec(END_CALL, fa, 1)
    s.rec(END_SCOPE, ib)
    s.rec(END_SCOPE, ia)
    s.rec(END_CALL, m, 1)
    s.rec(END_CALL, root, 1)
    S = ca.StatsRecorder.from_bytes(s.buf)
    self.assertFalse(S.truncated())
    roots = S.to_native()
    self.assertEqual(len(roots), 1)
    G = roots[0]
    self.assertEqual(G["name"], "G")
    m = G["children"][0]
    self.assertNotIn("children", m)
    [sb], [sa] = [it["children"] for it in m["iterations"]]
    self.assertEqual(sa["flag"], 1)
    self.assertEqual((sa["mem"], sb["mem"]), (0, 1))
    self.assertEqual(sa["stats"]["iterations"], {"iter": [0, 1], "obj": [1.5, 1.0], "mu": [0.1, 0.01]})
    self.assertEqual(sb["stats"]["iterations"], {"iter": [0], "obj": [2.5], "mu": [0.2]})
    self.assertEqual(sa["stats"]["warnings"], ["restoration"])
    self.assertEqual(sa["stats"]["iter_count"], 2)
    self.assertEqual(sb["stats"]["return_status"], "Solve_Succeeded")
    self.assertEqual(sa["id"], "p_f2:solver")
    # Iterations without calls; no sections were opened
    self.assertEqual(sa["iterations"], [{"children": []}] * 2)
    self.assertEqual(sorted(sa), ["flag", "id", "iterations", "mem", "name", "stats"])
    # Last call in stream order
    self.assertEqual(S.get_stat("solver", "return_status"), "Solve_Succeeded")
    with self.assertRaises(Exception): S.get_stat("solver", "iter_count")
    j = json.loads(S.to_json())
    self.assertEqual(j["version"], 1)
    self.assertEqual(j["calls"][0]["children"][0]["iterations"][1]["children"][0]["flag"], 1)
    U = S.deinterleave()
    self.assertEqual(U.to_native(), roots)
    with tempfile.TemporaryDirectory() as d:
      S.save(os.path.join(d, "s.cbor"))
      self.assertEqual(ca.StatsRecorder.load(os.path.join(d, "s.cbor")).to_native(), roots)
      # .json is for export_json only; distinct names for case-insensitive file systems
      for f in ["a.json", "b.JSON"]:
        with self.assertRaises(Exception): S.save(os.path.join(d, f))
        self.assertFalse(os.path.exists(os.path.join(d, f)))
        S.export_json(os.path.join(d, f))
        with open(os.path.join(d, f)) as h: self.assertEqual(h.read(), S.to_json())
        with self.assertRaises(Exception): ca.StatsRecorder.load(os.path.join(d, f))

  def test_decode_sections(self):
    s = _Stream()
    root = s.rec(BEGIN_CALL, None, "p_f0:solver", None)
    call = lambda parent, n: s.rec(END_CALL, s.rec(BEGIN_CALL, parent, "p_f1:" + n, None), 0)
    s.rec(DECLARE_FIELDS, root, ["mode", "obj", "trials"])
    pre = s.rec(BEGIN_SECTION, root, "pre")
    call(pre, "init")
    s.rec(END_SCOPE, pre)
    it = s.rec(BEGIN_ITERATION, root, 0)
    s.rec(SET_FIELD, it, 0, "regular")
    s.rec(SET_FIELD, it, 1, 3.0)
    s.rec(SET_FIELD, it, 2, 0)
    call(it, "a0")
    call(it, "b0")
    s.rec(END_SCOPE, it)
    it = s.rec(BEGIN_ITERATION, root, 1)
    s.rec(SET_FIELD, it, 0, "restoration")
    s.rec(SET_FIELD, it, 1, 2)
    s.rec(SET_FIELD, it, 2, 1)
    s.rec(SET, it, "warning", "tiny step")
    call(it, "a1")
    s.rec(END_SCOPE, it)
    # obj and trials missing; a field that was not declared
    it = s.rec(BEGIN_ITERATION, root, 2)
    s.rec(SET_FIELD, it, 0, "regular")
    s.rec(SET_FIELD, it, 3, True)
    s.rec(END_SCOPE, it)
    post = s.rec(BEGIN_SECTION, root, "post")
    call(post, "grad")
    s.rec(END_SCOPE, post)
    s.rec(SET, root, "iter_count", 2)
    s.rec(END_CALL, root, 0)
    S = ca.StatsRecorder.from_bytes(s.buf)
    [r] = S.to_native()
    self.assertEqual(sorted(r), ["flag", "id", "iterations", "mem", "name", "post", "pre", "stats"])
    for sec in [r["pre"], r["post"], r["iterations"][0], r["iterations"][2]]:
      self.assertEqual(sorted(sec), ["children"])
    self.assertEqual(names(r["pre"]), ["init"])
    self.assertEqual([names(it) for it in r["iterations"]], [["a0", "b0"], ["a1"], []])
    self.assertEqual(r["iterations"][1]["stats"], {"warning": "tiny step"})
    self.assertEqual(names(r["post"]), ["grad"])
    # Typed columns; missing: NaN if numeric, else None
    cols = r["stats"]["iterations"]
    self.assertEqual(sorted(cols), ["field3", "iter", "mode", "obj", "trials"])
    self.assertEqual(cols["iter"], [0, 1, 2])
    self.assertEqual(cols["mode"], ["regular", "restoration", "regular"])
    self.assertEqual(cols["obj"][:2], [3.0, 2.0])
    self.assertEqual(cols["trials"][:2], [0, 1])
    self.assertTrue(np.isnan(cols["obj"][2]) and np.isnan(cols["trials"][2]))
    self.assertEqual(cols["field3"], [None, None, True])
    self.assertEqual(r["stats"]["iter_count"], 2)
    # (NaN compares unequal to itself; JSON writes it as null)
    self.assertEqual(S.deinterleave().to_json(), S.to_json())
    self.assertEqual(json.loads(S.to_json())["calls"][0]["iterations"], r["iterations"])
    for n in s.bounds:
      [c] = ca.StatsRecorder.from_bytes(s.buf[:n]).to_native() or [r]
      self.assertTrue(is_prefix(c, r))
    s = _Stream()
    root = s.rec(BEGIN_CALL, None, "p_f0:F", None)
    sec = s.rec(BEGIN_SECTION, root, "setup")
    call(sec, "f")
    s.rec(END_SCOPE, sec)
    s.rec(END_CALL, root, 0)
    [r] = ca.StatsRecorder.from_bytes(s.buf).to_native()
    self.assertEqual(sorted(r), ["flag", "id", "mem", "name", "setup", "stats"])
    self.assertEqual(names(r["setup"]), ["f"])
    self.assertEqual(names(r["setup"]["children"][0]), [])

  def test_decode_truncated(self):
    s = _Stream()
    r = s.rec(BEGIN_CALL, None, "p_f0:F", None)
    for i in range(5):
      s.rec(END_CALL, s.rec(BEGIN_CALL, r, "p_f1:f", None), 0)
    s.rec(SET, r, "x", 1.25)
    s.rec(END_CALL, r, 0)
    full = ca.StatsRecorder.from_bytes(s.buf).to_native()[0]
    for n in range(len(s.buf) + 1):
      # Cut with 0xff and stale bytes, or plain mid-record end
      tails = [b"", b"\xff", b"\xff\x00\x81\x02"] if n in s.bounds else [b""]
      for tail in tails:
        S = ca.StatsRecorder.from_bytes(s.buf[:n] + tail)
        roots = S.to_native()
        self.assertEqual(S.truncated(), bool(tail) or n not in s.bounds)
        self.assertTrue(all(is_prefix(c, full) for c in roots))
        if n < len(s.buf): self.assertFalse(roots and roots[0]["flag"] is not None)

  # --- virtual machine ---

  def test_vm_tree(self):
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y")
    g = ca.Function("g", [y], [f(y) + f(2 * y)])
    S = ca.StatsRecorder()
    g(S, 0.3)
    roots = S.to_native()
    self.assertEqual(shape(roots[0]), ("g", 0, (("children", ("f", 0, ())),) * 2))
    # Equal names, distinct ids
    self.assertEqual(len({c["id"] for c in find(roots, "f")}), 1)
    self.assertNotEqual(find(roots, "f")[0]["id"], roots[0]["id"])
    g(S, 0.3)
    self.assertEqual(len(S.to_native()), 1)
    S.clear()
    self.assertEqual((S.nbytes(), S.to_native(), S.truncated()), (0, [], False))
    g(S, 0.3)
    self.assertEqual(shape(S.to_native()[0]), shape(roots[0]))
    nbytes = S.nbytes()
    g(0.3)
    self.assertEqual(S.nbytes(), nbytes)
    # Both call forms record
    S2 = ca.StatsRecorder()
    g.call(S2, {"i0": 0.3})
    self.assertEqual(shape(S2.to_native()[0]), shape(roots[0]))

  def test_vm_gather_stats(self):
    # gather_stats off: neither the call nor its callees
    x = ca.SX.sym("x")
    y = ca.MX.sym("y")
    f = ca.Function("f", [x], [ca.sin(x)])
    for opts, expected in [({}, 2), ({"gather_stats": True}, 2), ({"gather_stats": False}, 0)]:
      h = ca.Function("h", [y], [f(y)], opts)
      g = ca.Function("g", [y], [h(2 * y) + h(3 * y)])
      # Also after deserialization
      for G in [g, ca.Function.deserialize(g.serialize())]:
        S = ca.StatsRecorder()
        self.checkarray(G(S, 0.3), np.sin(0.6) + np.sin(0.9))
        [r] = S.to_native()
        self.assertEqual(len(find([r], "h")), expected)
        self.assertEqual(len(find([r], "f")), expected)
    # At the top, nothing is recorded
    S = ca.StatsRecorder()
    ca.Function("h", [y], [f(2 * y)], {"gather_stats": False})(S, 0.3)
    self.assertEqual((S.nbytes(), S.to_native()), (0, []))

  def test_vm_capacity(self):
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    n = 10000
    y = ca.MX.sym("y", 1, n)
    g = ca.Function("g", [y], [f.map(n)(y)])
    L = ca.StatsRecorder({"size": 1 << 20})
    g(L, np.arange(float(n)))
    self.assertTrue(L.nbytes() > 1 << 17)
    self.assertFalse(L.truncated())
    self.assertEqual(len(L.to_native()[0]["children"][0]["iterations"]), n)
    # Truncated
    M = ca.StatsRecorder({"size": 4096})
    g(M, np.arange(float(n)))
    self.assertTrue(M.truncated())
    self.assertEqual(M.nbytes(), L.nbytes())
    self.assertTrue(is_prefix(M.to_native()[0], L.to_native()[0]))

  def test_vm_map(self):
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y", 1, 4)
    for par in ["serial", "openmp", "thread"]:
      F = f.map(4, par)
      g = ca.Function("g", [y], [F(y)])
      S = ca.StatsRecorder()
      g(S, np.arange(4.))
      m = S.to_native()[0]["children"][0]
      self.assertEqual(m["name"], F.name())
      self.assertEqual(sorted(m), ["flag", "id", "iterations", "mem", "name", "stats"])
      self.assertEqual([[shape(c) for c in it["children"]] for it in m["iterations"]],
                       [[("f", 0, ())]] * 4)
      self.assertEqual(sorted(m["iterations"][0]), ["children"])

  def test_vm_dump(self):
    x = ca.SX.sym("x")
    with tempfile.TemporaryDirectory() as d:
      f = ca.Function("f", [x], [ca.sin(x)], {"dump_in": True, "dump_out": True, "dump": True,
                                              "dump_trace": True, "dump_dir": d})
      y = ca.MX.sym("y")
      g = ca.Function("g", [y], [f(y) + f(2 * y)])
      S = ca.StatsRecorder()
      g(S, 0.3)
      a, b = find(S.to_native(), "f")
      self.assertEqual(sorted(a["stats"]), ["dump", "dump_in", "dump_out", "dump_trace"])
      self.assertEqual(a["stats"]["dump"], os.path.join(d, "f.casadi"))
      for c, i in [(a, 0), (b, 1)]:
        for k, e in [("dump_in", "in.txt"), ("dump_out", "out.txt"), ("dump_trace", "trace.jsonl")]:
          self.assertEqual(c["stats"][k], os.path.join(d, "f.%06d.%s" % (i, e)))
          self.assertTrue(os.path.exists(c["stats"][k]))
      self.assertNotIn("dump", b["stats"])
      self.assertTrue(os.path.exists(a["stats"]["dump"]))
    self.assertEqual(find(S.to_native(), "g")[0]["stats"], {})

  @requires_nlpsol("ipopt")
  def test_vm_ipopt(self):
    R = _Rosenbrock()
    S = ca.StatsRecorder()
    R.solver(S, p=100, **R.args)
    ref = R.solver.stats()
    s = S.to_native()[0]
    for k in ["return_status", "unified_return_status", "success", "iter_count"]:
      self.assertEqual(s["stats"][k], ref[k])
    for k, v in ref["iterations"].items():
      self.assertEqual(s["stats"]["iterations"][k], list(v))
    self.assertEqual(S.get_stat("solver", "return_status"), ref["return_status"])
    self.assertEqual(S.get_stat("solver", "iter_count"), ref["iter_count"])
    self.assertEqual((s["stats"]["return_status_enum"], s["stats"]["unified_return_status_enum"]),
                     (0, 0))
    # Maximum_Iterations_Exceeded (-1), SOLVER_RET_LIMITED (2)
    x = ca.MX.sym("x", 2)
    limited = ca.nlpsol("limited", "ipopt", {"x": x, "f": ca.sumsqr(x - 1) + x[0]**4},
      {"print_time": False, "ipopt.print_level": 0, "ipopt.sb": "yes", "ipopt.max_iter": 1})
    S = ca.StatsRecorder()
    limited(S, x0=[5, 5])
    ls = S.to_native()[0]["stats"]
    self.assertEqual((ls["return_status"], ls["return_status_enum"]),
                     ("Maximum_Iterations_Exceeded", -1))
    self.assertEqual((ls["unified_return_status"], ls["unified_return_status_enum"]),
                     ("SOLVER_RET_LIMITED", 2))
    self.assertEqual(len(s["iterations"]), ref["iter_count"] + 1)
    self.assertTrue({"nlp_f", "nlp_g", "nlp_grad_f", "nlp_jac_g"} <= set(names(s["pre"])))
    self.assertNotIn("nlp_hess_l", names(s["pre"]))
    for it in s["iterations"][:-1]:
      self.assertEqual(names(it)[0], "nlp_hess_l")
      self.assertIn("nlp_f", names(it))
    self.assertEqual(s["iterations"][-1], {"children": []})
    self.assertEqual(names(s["post"]), ["nlp_grad"])
    cols = s["stats"]["iterations"]
    self.assertEqual(cols["iter"], list(range(ref["iter_count"] + 1)))
    self.assertEqual(set(cols["alg_mod"]), {"regular"})
    self.assertTrue(all(isinstance(e, int) for e in cols["ls_trials"]))
    self.assertEqual(ref["n_call_nlp_f"], names(calls(s)).count("nlp_f"))
    self.assertEqual(ref["n_call_nlp_hess_l"], names(calls(s)).count("nlp_hess_l"))

  @requires_nlpsol("sqpmethod")
  def test_vm_nlpsol_flat(self):
    # Plugin without iteration recording: oracle calls are plain children
    x = ca.MX.sym("x", 2)
    solver = ca.nlpsol("solver", "sqpmethod", {"x": x, "f": ca.sumsqr(x - 1)},
      {"print_time": False, "print_iteration": False, "print_header": False,
       "print_status": False, "qpsol": "qrqp",
       "qpsol_options": {"print_iter": False, "print_header": False, "print_info": False}})
    S = ca.StatsRecorder()
    solver(S, x0=[0, 0])
    s = S.to_native()[0]
    self.assertEqual(sorted(s), ["children", "flag", "id", "mem", "name", "stats"])
    self.assertIn("nlp_hess_l", names(s["children"]))

  @requires_nlpsol("ipopt")
  def test_vm_ipopt_callback_step(self):
    # iteration_callback_step thins the callback, not the iterations recorded
    class Counter(ca.Callback):
      def __init__(self, nx, ng):
        ca.Callback.__init__(self)
        self.nx, self.ng, self.iters = nx, ng, []
        self.construct("counter", {})
      def get_n_in(self): return ca.nlpsol_n_out()
      def get_n_out(self): return 1
      def get_name_in(self, i): return ca.nlpsol_out(i)
      def get_name_out(self, i): return "ret"
      def get_sparsity_in(self, i):
        n = ca.nlpsol_out(i)
        if n in ("f", "lam_p", "p"): return ca.Sparsity.dense(1 if n == "f" else 0)
        return ca.Sparsity.dense(self.nx if n in ("x", "lam_x") else self.ng)
      def eval(self, arg):
        self.iters.append(float(arg[0][0]))
        return [0]
    R = _Rosenbrock()
    V = ca.StatsRecorder()
    R.solver(V, p=100, **R.args)
    ref = V.to_native()[0]
    nlp = R.solver.oracle()
    for step in [1, 2, 3]:
      cb = Counter(2, 1)
      solver = ca.nlpsol("solver", "ipopt", nlp,
        {"print_time": False, "ipopt.print_level": 0, "ipopt.sb": "yes",
         "iteration_callback": cb, "iteration_callback_step": step})
      S = ca.StatsRecorder()
      solver(S, p=100, **R.args)
      s = S.to_native()[0]
      n = solver.stats()["iter_count"]
      self.assertEqual(len(cb.iters), len(range(0, n + 1, step)))
      self.assertEqual(len(solver.stats()["iterations"]["obj"]), n + 1)
      self.assertEqual(s["stats"]["iterations"], ref["stats"]["iterations"])
      self.assertEqual(shape(s), shape(ref))

  @requires_nlpsol("ipopt")
  def test_vm_deinterleave(self):
    R = _Rosenbrock()
    P = ca.MX.sym("P", 1, 4)
    H = ca.Function("H", [P], [R.solver.map(4, "thread")(
      x0=ca.repmat(ca.DM([-1, 1]), 1, 4), p=P, lbg=-ca.inf, ubg=1)["x"]])
    S = ca.StatsRecorder({"size": 1e6})
    H(S, [1., 10., 100., 1000.])
    U = S.deinterleave()
    self.assertEqual([strip_ids(r) for r in U.to_native()], [strip_ids(r) for r in S.to_native()])
    m = S.to_native()[0]["children"][0]
    self.assertEqual(len(m["iterations"]), 4)
    for s in [s for it in m["iterations"] for s in it["children"]]:
      self.assertEqual(s["name"], "solver")
      self.assertEqual(len(s["iterations"]), s["stats"]["iter_count"] + 1)
      self.assertEqual(names(s["post"]), ["nlp_grad"])

  # --- generated code ---

  def test_cli(self):
    # casadi-cli stats
    cli = os.path.join(ca.GlobalOptions.getCasadiPath(),
                       "casadi-cli" + (".exe" if os.name == "nt" else ""))
    if not os.path.exists(cli): self.skipTest("casadi-cli not found")
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y", 1, 3)
    g = ca.Function("g", [y], [f.map(3, "thread")(y)])
    S = ca.StatsRecorder()
    g(S, [1, 2, 3])
    run = lambda *a: subprocess.run([cli, "stats"] + list(a), capture_output=True, text=True)
    with tempfile.TemporaryDirectory() as d:
      fn = os.path.join(d, "g.cbor")
      S.save(fn)
      r = run(fn)
      self.assertEqual(r.returncode, 0)
      self.assertEqual(r.stdout.strip().splitlines()[1:], S.str(True).strip().splitlines()[1:])
      self.assertEqual([l.strip() for l in r.stdout.splitlines()[1:]],
                       ["g", "threadmap3_f", "iteration 0:", "f", "iteration 1:", "f",
                        "iteration 2:", "f"])
      r = run(fn, "--json")
      self.assertEqual((r.returncode, r.stdout.strip()), (0, S.to_json()))
      out = os.path.join(d, "g.json")
      self.assertEqual(run(fn, "--json", out).returncode, 0)
      with open(out) as j: self.assertEqual(j.read(), S.to_json())
      # Truncated stream: exit code 2
      if ca.has_nlpsol("ipopt"):
        R = _Rosenbrock()
        T = ca.StatsRecorder()
        R.solver(T, **R.args)
        T.save(fn)
        r = run(fn, "--get", "solver", "return_status")
        self.assertEqual((r.returncode, r.stdout.strip()), (0, "Solve_Succeeded"))
        self.assertIn("return_status: Solve_Succeeded", run(fn).stdout)
      T = ca.StatsRecorder({"size": 40})
      g(T, [1, 2, 3])
      T.save(fn)
      self.assertEqual(run(fn).returncode, 2)
      r = run(fn, "--bogus")
      self.assertEqual(r.returncode, 1)
      self.assertIn("--json", r.stderr)

  @slow()
  def test_codegen_tree(self):
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y")
    g = ca.Function("g", [y], [f(y) + f(2 * y)])
    lib = self.build(g, [0.3])
    flag, outs, S, needed = self.call(lib, g, [0.3])
    self.assertEqual(flag, 0)
    self.checkarray(outs[0], g(0.3))
    V = ca.StatsRecorder()
    g(V, 0.3)
    root = S.to_native()[0]
    self.assertEqual(shape(root), shape(V.to_native()[0]))
    self.assertTrue(root["id"].split(":")[0].endswith("_f0"))
    # Owner reserve, mimicking the default; commit sees each record complete, once
    stream = self.stream
    asked = {}
    committed = []
    def reserve(s, n, pos):
      st = ctypes.cast(s, ctypes.POINTER(CStats)).contents
      pos[0] = st.needed
      st.needed += n
      asked[pos[0]] = n
      return ctypes.addressof(st.p.contents) + pos[0] if pos[0] + n <= st.cap else None
    def commit(s, pos):
      st = ctypes.cast(s, ctypes.POINTER(CStats)).contents
      committed.append((pos, ctypes.string_at(ctypes.addressof(st.p.contents) + pos, asked[pos])))
    flag, outs2, S, needed = self.call(lib, g, [0.3], reserve=CReserve(reserve),
                                       commit=CCommit(commit))
    self.assertEqual(flag, 0)
    self.assertEqual(self.stream, stream)
    self.assertEqual(sum(asked.values()), len(stream))
    self.assertEqual(sorted(pos for pos, _ in committed), sorted(asked))
    for pos, rec in committed: self.assertEqual(rec, stream[pos:pos + asked[pos]])
    # Commit alone, with the default reserve
    starts = []
    flag, outs2, S, needed = self.call(lib, g, [0.3], commit=CCommit(lambda s, pos: starts.append(pos)))
    self.assertEqual(flag, 0)
    self.assertEqual(self.stream, stream)
    self.assertEqual(starts, sorted(asked))
    flag, outs2, S, needed = self.call(lib, g, [0.3], raw=True)
    self.assertEqual(flag, 0)
    self.assertEqual(self.stream, stream)
    flag, outs2, S, needed = self.call(lib, g, [0.3], with_stats=False)
    self.assertEqual(needed, 0)
    self.checkarray(outs2[0], outs[0])

  @slow()
  def test_codegen_unrolled(self):
    # <name>_unrolled_with_stats records the same stream as <name>_with_stats
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y")
    z = ca.MX.sym("z", 2)
    g = ca.Function("g", [y, z], [f(y) + f(2 * y), z * y])
    lib = self.build(g, [0.3, [1, 2]], opts={"unroll_args": True})
    flag, outs, S, needed = self.call(lib, g, [0.3, [1, 2]])
    stream = self.stream
    flag, outs2, S, needed = self.call(lib, g, [0.3, [1, 2]], unrolled=True)
    self.assertEqual(flag, 0)
    self.assertEqual(self.stream, stream)
    for a, b in zip(outs2, g(0.3, [1, 2])): self.checkarray(a, b)
    # Without stats, no unrolled stats entry point
    cg = ca.CodeGenerator("a", {"unroll_args": True}); cg.add(g)
    self.assertIn("g_unrolled(", cg.dump())
    self.assertNotIn("_unrolled_with_stats", cg.dump())

  @slow()
  def test_codegen_foreign_reserve(self):
    # Without function pointers, an owner reserve cannot be called: records dropped
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y")
    g = ca.Function("g", [y], [f(y) + f(2 * y)])
    lib = self.build(g, [0.3], opts={"allow_function_pointers": False})
    asked = []
    def reserve(s, n, pos):
      asked.append(n)
      return None
    flag, outs, S, needed = self.call(lib, g, [0.3], reserve=CReserve(reserve))
    self.assertEqual(flag, 0)
    self.checkarray(outs[0], g(0.3))
    self.assertEqual((asked, needed, self.stream), ([], 0, b""))
    # Likewise an owner commit
    committed = []
    flag, outs, S, needed = self.call(lib, g, [0.3], commit=CCommit(lambda s, pos: committed.append(pos)))
    self.assertEqual(flag, 0)
    self.assertEqual((committed, needed, self.stream), ([], 0, b""))

  @slow()
  def test_codegen_reserve_exported(self):
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    for fp in [True, False]:
      lib = self.build(f, [0.3], opts={"allow_function_pointers": fp})
      reserve = getattr(lib, getattr(lib, "_prefix") + "_stats_reserve")
      reserve.restype = ctypes.c_void_p
      buf = (ctypes.c_ubyte * 16)()
      sink = CStats(ctypes.cast(buf, ctypes.POINTER(ctypes.c_ubyte)), 16, 0, CReserve(), CCommit(), None)
      pos = ctypes.c_longlong()
      base = ctypes.addressof(buf)
      self.assertEqual(reserve(ctypes.byref(sink), ctypes.c_longlong(10), ctypes.byref(pos)), base)
      self.assertEqual((pos.value, sink.needed), (0, 10))
      self.assertEqual(reserve(ctypes.byref(sink), ctypes.c_longlong(4), ctypes.byref(pos)), base + 10)
      self.assertEqual((pos.value, sink.needed), (10, 14))
      # Past the end: null, 0xff at the cut, needed keeps counting
      self.assertIsNone(reserve(ctypes.byref(sink), ctypes.c_longlong(5), ctypes.byref(pos)))
      self.assertEqual((pos.value, sink.needed, buf[14]), (14, 19, 0xff))
      self.assertIsNone(reserve(ctypes.byref(sink), ctypes.c_longlong(1), ctypes.byref(pos)))
      self.assertEqual((pos.value, sink.needed), (19, 20))
      nbytes = getattr(lib, getattr(lib, "_prefix") + "_stats_nbytes")
      nbytes.restype = ctypes.c_longlong
      truncated = getattr(lib, getattr(lib, "_prefix") + "_stats_truncated")
      data = getattr(lib, getattr(lib, "_prefix") + "_stats_data")
      data.restype = ctypes.c_void_p
      n = ctypes.c_longlong()
      self.assertEqual((nbytes(ctypes.byref(sink)), truncated(ctypes.byref(sink))), (20, 1))
      self.assertEqual((data(ctypes.byref(sink), ctypes.byref(n)), n.value), (base, 16))
      getattr(lib, getattr(lib, "_prefix") + "_stats_clear")(ctypes.byref(sink))
      self.assertEqual(sink.needed, 0)
      self.assertEqual((nbytes(ctypes.byref(sink)), truncated(ctypes.byref(sink))), (0, 0))
      self.assertEqual(reserve(ctypes.byref(sink), ctypes.c_longlong(3), ctypes.byref(pos)), base)
      make = getattr(lib, getattr(lib, "_prefix") + "_stats_make_sink")
      make.restype = CStats
      made = make(ctypes.cast(buf, ctypes.POINTER(ctypes.c_ubyte)), ctypes.c_longlong(16))
      self.assertEqual((ctypes.addressof(made.p.contents), made.cap, made.needed), (base, 16, 0))
      self.assertEqual(bool(made.reserve), fp)
      self.assertFalse(made.commit)
      self.assertEqual(reserve(ctypes.byref(made), ctypes.c_longlong(6), ctypes.byref(pos)), base)
      self.assertEqual((pos.value, made.needed), (0, 6))

  @slow()
  def test_codegen_gather_stats(self):
    # At call sites and at the entry point
    x = ca.SX.sym("x")
    y = ca.MX.sym("y")
    f = ca.Function("f", [x], [ca.sin(x)])
    h = ca.Function("h", [y], [f(2 * y)], {"gather_stats": False})
    g = ca.Function("g", [y], [h(2 * y) + h(3 * y) + f(4 * y)])
    lib = self.build(g, [0.3])
    flag, outs, S, needed = self.call(lib, g, [0.3])
    self.assertEqual(flag, 0)
    self.checkarray(outs[0], g(0.3))
    V = ca.StatsRecorder()
    g(V, 0.3)
    self.assertEqual(shape(S.to_native()[0]), shape(V.to_native()[0]))
    self.assertEqual(names(S.to_native()[0]["children"]), ["f"])
    lib = self.build(h, [0.3])
    flag, outs, S, needed = self.call(lib, h, [0.3])
    self.assertEqual((flag, needed), (0, 0))
    self.checkarray(outs[0], np.sin(0.6))

  @slow()
  def test_codegen_map(self):
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y", 1, 4)
    for par, fp in itertools.product(["serial", "openmp", "thread"], [True, False]):
      g = ca.Function("g", [y], [f.map(4, par)(y) * 2])
      flags, defs = parallel_options(par)
      lib = self.build(g, [np.arange(4.)],
                       opts={"thread_safe": par != "serial", "allow_function_pointers": fp},
                       std="c99", extra_options=flags, definitions=defs, main=True,
                       valgrind=True, helgrind=par == "thread")
      flag, outs, S, needed = self.call(lib, g, [np.arange(4.)])
      self.assertEqual(flag, 0)
      self.checkarray(outs[0], 2 * np.sin(np.arange(4.)))
      V = ca.StatsRecorder()
      g(V, np.arange(4.))
      self.assertEqual(shape(S.to_native()[0]), shape(V.to_native()[0]))

  @slow()
  def test_codegen_map_mutex(self):
    # Mutex fallback, with a library called from the iterations; C11: mutex set up in incref
    if os.name == "nt": self.skipTest("the fallback is reached by hiding the GNU atomic")
    x = ca.SX.sym("x")
    y = ca.MX.sym("y", 1, 16)

    def compile_without_atomic(F, opts, definitions, extralibs=None):
      # Hide the GNU atomic
      cg = ca.CodeGenerator(F.name(), dict(opts, stats=True, with_header=True))
      cg.add(F)
      cg.generate()
      with open(F.name() + ".c") as src:
        code = src.read().replace("#elif defined(__GNUC__) || defined(__clang__)", "#elif 0")
      self.assertIn("CASADI_MUTEX_LOCK(&casadi_stats_mutex);", code)
      with open(F.name() + ".c", "w") as src: src.write(code)
      return self.compile_external(F.name(), F.name() + ".c", std="c11",
                                   extra_options=["-pthread"], definitions=definitions,
                                   extralibs=extralibs or "")

    # C11 threads.h: Linux only
    for fp, threads in itertools.product([True, False],
                                         ["POSIX"] + (["C11"] if sys.platform.startswith("linux") else [])):
      opts = {"thread_safe": True, "allow_function_pointers": fp}
      defs = ["CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_" + threads]
      uid = hashlib.md5((str(opts) + threads + str(time.time())).encode()).hexdigest()
      f = ca.Function("fmutex_" + uid, [x], [ca.sin(x)])
      F, fso = compile_without_atomic(f, opts, defs)
      # Without function pointers, a called library would take its own mutex: records may collide
      for inner in [f, F] if fp else [f]:
        g = ca.Function("gmutex_%s_%d" % (uid, inner is F), [y], [inner.map(16, "thread")(y) * 2])
        _, gso = compile_without_atomic(g, opts, defs, [fso] if inner is F else None)
        lib = cdll(gso, g.name())
        for rep in range(5):
          flag, outs, S, needed = self.call(lib, g, [np.arange(16.)])
          self.assertEqual(flag, 0)
          self.checkarray(outs[0], 2 * np.sin(np.arange(16.)))
          self.assertEqual(needed, len(self.stream))
          self.assertFalse(S.truncated())
          if inner is f:
            V = ca.StatsRecorder()
            g(V, np.arange(16.))
            self.assertEqual(shape(S.to_native()[0]), shape(V.to_native()[0]))
          else:
            [m] = [c for c in walk(S.to_native()[0]) if "iterations" in c]
            self.assertEqual(len(m["iterations"]), 16)
            for it in m["iterations"]:
              [c] = it["children"]
              [d] = c["children"]
              self.assertEqual((c["name"], d["name"]), (f.name(), f.name()))
              self.assertTrue(d["id"].startswith(f.name()))

  @slow()
  def test_codegen_vm_parity(self):
    # Same call tree from VM and generated code: switch, FD, mapsum, SX call node
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    h = ca.Function("h", [x], [ca.cos(x)])
    c = ca.MX.sym("c")
    y = ca.MX.sym("y")
    sw = ca.Function.conditional("sw", [f], h)
    F = ca.Function("F", [x], [ca.sin(x)], {"enable_fd": True, "enable_forward": False,
                                           "enable_reverse": False, "enable_jacobian": False})
    fd = F.forward(1)
    v = ca.MX.sym("v", 1, 3)
    z = ca.SX.sym("z")
    k = ca.Function("k", [z], [f.call([z], False, True)[0] * 2])
    cases = [(ca.Function("g", [c, y], [sw(c, y)]), [[0, 0.3], [1, 0.3]]),
             (ca.Function("g", [y], [fd(y, F(y), 1)]), [[0.3]]),
             (ca.Function("g", [v], [f.map(3, [False], [True])(v)]), [[ca.DM([[0.1, 0.2, 0.3]])]]),
             (ca.Function("g", [y], [k(y)]), [[0.3]])]
    for g, inputs in cases:
      # Finite differences need C99 (NAN, fmax)
      lib = self.build(g, inputs[0], std="c99")
      for inp in inputs:
        flag, outs, S, needed = self.call(lib, g, inp)
        self.assertEqual(flag, 0)
        self.checkarray(outs[0], g(*inp))
        V = ca.StatsRecorder()
        g(V, *inp)
        self.assertEqual(shape(S.to_native()[0]), shape(V.to_native()[0]))
        self.assertTrue(len(list(walk(V.to_native()[0]))) > 2)

  @requires_nlpsol("ipopt")
  @slow()
  def test_codegen_ipopt(self):
    R = _Rosenbrock()
    S = R.solver
    inputs = [[-1, 1], 100, -ca.inf, ca.inf, -ca.inf, 1, 0, 0]
    lib = self.build(S, inputs, std="c99", extralibs=["ipopt"], digits=8)
    flag, outs, C, needed = self.call(lib, S, inputs)
    self.assertEqual(flag, 0)
    s = C.to_native()[0]
    V = ca.StatsRecorder()
    S(V, p=100, **R.args)
    ref = S.stats()
    v = V.to_native()[0]
    for k in ["return_status", "unified_return_status", "success", "iter_count"]:
      self.assertEqual(s["stats"][k], ref[k])
    # Same iterates, bit for bit
    for k, e in ref["iterations"].items():
      self.assertEqual(s["stats"]["iterations"][k], list(e))
    self.assertEqual(shape(s), shape(v))
    self.assertEqual(len(s["iterations"]), ref["iter_count"] + 1)
    self.assertEqual(names(s["post"]), ["nlp_grad"])
    # Byte-identical to generated code
    self.assertEqual(s["stats"], v["stats"])
    self.assertEqual(len(s["stats"]["iterations"]["iter"]), ref["iter_count"] + 1)

  @requires_nlpsol("ipopt")
  @slow()
  def test_codegen_ipopt_truncated(self):
    R = _Rosenbrock()
    S = R.solver
    inputs = [[-1, 1], 100, -ca.inf, ca.inf, -ca.inf, 1, 0, 0]
    lib = self.build(S, inputs, std="c99", extralibs=["ipopt"], digits=8)
    _, xref, full, n = self.call(lib, S, inputs)
    ref = full.to_native()[0]
    for cap in list(range(0, 64)) + list(range(64, n + 2, 37)) + [n - 1, n, n + 1]:
      flag, outs, C, needed = self.call(lib, S, inputs, cap=cap)
      self.assertEqual(flag, 0)
      self.checkarray(outs[0], xref[0])
      self.assertEqual(needed, n)
      if cap >= n: self.assertFalse(C.truncated())
      self.assertTrue(all(is_prefix(r, ref) for r in C.to_native()))

  @slow()
  def test_codegen_get_stat_types(self):
    x = ca.SX.sym("x")
    lib = self.build(ca.Function("f", [x], [ca.sin(x)]), [0.3])
    s = _Stream()
    root = s.rec(BEGIN_CALL, None, "p_f0:G", None)
    a = s.rec(BEGIN_CALL, root, "p_f1:solver", 0)
    s.rec(SET, a, "old", 7)
    s.rec(END_CALL, a, 0)
    b = s.rec(BEGIN_CALL, root, "p_f1:solver", 0)
    s.rec(DECLARE_FIELDS, b, ["n"])
    pre = s.rec(BEGIN_SECTION, b, "pre")
    s.rec(SET, pre, "n", 97)
    s.rec(END_SCOPE, pre)
    it = s.rec(BEGIN_ITERATION, b, 0)
    s.rec(SET, it, "n", 98)
    s.rec(SET_FIELD, it, 0, 99)
    s.rec(END_SCOPE, it)
    for k, v in [("n", 3), ("neg", -5), ("big", 2**40), ("yes", True), ("no", False), ("t", 0.25),
                 ("status", "Solve_Succeeded"), ("names", ["obj", "mu"]), ("none", None)]:
      s.rec(SET, b, k, v)
    s.rec(SET, b, "n", 4)
    s.rec(APPEND, b, "row", [1.0, 2.0])
    s.rec(SET, root, "outer", 9)
    s.rec(END_CALL, b, 0)
    s.rec(END_CALL, root, 0)
    buf = s.buf
    M = "missing"
    for key, i, r, t in [("n", 4, 4.0, M), ("neg", -5, -5.0, M), ("big", 2**40, 2.0**40, M),
                         ("yes", 1, 1.0, M), ("no", 0, 0.0, M), ("t", M, 0.25, M),
                         ("status", M, M, "Solve_Succeeded"), ("names", M, M, M), ("none", M, M, M),
                         ("row", M, M, M), ("old", M, M, M)]:
      self.assertEqual(get_stat(lib, buf, "solver", key, "int"), i)
      self.assertEqual(get_stat(lib, buf, "solver", key, "real"), r)
      self.assertEqual(get_stat(lib, buf, "solver", key, "text"), t)
    # Truncated, null-terminated; cap 0 leaves dst alone
    self.assertEqual(get_stat(lib, buf, "solver", "status", "text", cap=6), "Solve")
    self.assertEqual(get_stat(lib, buf, "solver", "status", "text", cap=1), "")
    self.assertEqual(get_stat(lib, buf, "solver", "status", "text", cap=0), "#")
    self.assertEqual(get_stat(lib, buf, "G", "outer"), 9)
    self.assertEqual(get_stat(lib, buf, "olver", "n"), "missing")
    self.assertEqual(get_stat(lib, buf, "solver", "outer"), "missing")
    for cut in range(len(buf)):
      self.assertIn(get_stat(lib, buf[:cut], "solver", "n"), ["missing", 3, 4])
      self.assertIn(get_stat(lib, buf[:cut] + b"\xff" + buf[cut:], "solver", "n"), ["missing", 3, 4])
      self.assertIn(get_stat(lib, buf[:cut], "solver", "t", "real"), ["missing", 0.25])
      self.assertIn(get_stat(lib, buf[:cut], "solver", "status", "text"),
                    ["missing", "Solve_Succeeded"])

  @requires_nlpsol("ipopt")
  @slow()
  def test_codegen_get_stat(self):
    R = _Rosenbrock()
    S = R.solver
    inputs = [[-1, 1], 100, -ca.inf, ca.inf, -ca.inf, 1, 0, 0]
    lib = self.build(S, inputs, std="c99", extralibs=["ipopt"], digits=8)
    _, _, C, n = self.call(lib, S, inputs)
    stream = self.stream
    V = ca.StatsRecorder()
    S(V, p=100, **R.args)
    with tempfile.TemporaryDirectory() as d:
      V.save(os.path.join(d, "vm.cbor"))
      with open(os.path.join(d, "vm.cbor"), "rb") as f: vm = f.read()
    for s in [stream, vm]:
      self.assertEqual(get_stat(lib, s, "solver", "iter_count"), C.get_stat("solver", "iter_count"))
      self.assertEqual(get_stat(lib, s, "solver", "iter_count", "real"), C.get_stat("solver", "iter_count"))
      self.assertEqual(get_stat(lib, s, "solver", "success"), 1)
      self.assertEqual(get_stat(lib, s, "solver", "return_status_enum"), 0)
      self.assertEqual(get_stat(lib, s, "solver", "unified_return_status_enum"), 0)
      self.assertEqual(get_stat(lib, s, "solver", "return_status"), "missing")
      for k in ["return_status", "unified_return_status"]:
        self.assertEqual(get_stat(lib, s, "solver", k, "text"), C.get_stat("solver", k))
      self.assertEqual(get_stat(lib, s, "nlp_f", "success"), "missing")
    for cap in [0, 10, n // 2]:
      self.call(lib, S, inputs, cap=cap)
      self.assertEqual(get_stat(lib, self.stream, "solver", "success"), "missing")

  @requires_nlpsol("ipopt")
  @slow()
  def test_codegen_ipopt_threads(self):
    R = _Rosenbrock()
    P = ca.MX.sym("P", 1, 4)
    ps = [1., 10., 100., 1000.]
    for par, fp in itertools.product(["thread", "openmp"], [True, False]):
      H = ca.Function("H", [P], [R.solver.map(4, par)(
        x0=ca.repmat(ca.DM([-1, 1]), 1, 4), p=P, lbg=-ca.inf, ubg=1)["x"]])
      flags, defs = parallel_options(par)
      lib = self.build(H, [ps], opts={"thread_safe": True, "allow_function_pointers": fp},
                       std="c99", extra_options=flags, definitions=defs + ["CASADI_MAX_NUM_THREADS=4"],
                       extralibs=["ipopt"], digits=8)
      for rep in range(3):
        flag, outs, C, needed = self.call(lib, H, [ps], cap=1 << 20)
        self.assertEqual(flag, 0)
        self.assertFalse(C.truncated())
        roots = C.to_native()
        [m] = [c for r in roots for c in walk(r) if "iterations" in c and "pre" not in c]
        solves = [s for it in m["iterations"] for s in it["children"]]
        self.assertEqual(names(solves), ["solver"] * 4)
        # Memories may be reused with fewer threads than iterations
        self.assertTrue(all(0 <= s["mem"] < 4 for s in solves))
        for i, s in enumerate(solves):
          R.solver(p=ps[i], **R.args)
          ref = R.solver.stats()
          self.assertEqual(s["stats"]["iter_count"], ref["iter_count"])
          self.assertEqual(s["stats"]["iterations"]["obj"], list(ref["iterations"]["obj"]))
        self.assertEqual(C.deinterleave().to_native(), roots)

  @slow()
  def test_codegen_external(self):
    # Call tree continues into an external() generated library
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y")
    for lib_opts in [{}, {"stats": True}, {"stats": True, "allow_function_pointers": False}]:
      ret = self.check_codegen(f, inputs=[0.3], opts=dict(lib_opts))
      so = ret["libname"]
      name = os.path.splitext(os.path.basename(so))[0]
      lib_fp = lib_opts.get("stats", False) and lib_opts.get("allow_function_pointers", True)
      lib_stats = lib_opts.get("stats", False)
      F = ret["F"]
      g = ca.Function("g", [y], [F(y) + F(2 * y)])

      def check(roots, continued, id_prefix, fname="f"):
        g_ = roots[0]
        self.assertEqual([c["name"] for c in g_["children"]], [fname, fname])
        for c in g_["children"]:
          if continued:
            self.assertEqual(shape(c), (fname, 0, (("children", (fname, 0, ())),)))
            self.assertTrue(c["children"][0]["id"].startswith(id_prefix))
          else:
            self.assertEqual(c["children"], [])

      for fp in [True, False]:
        lib = self.build(g, [0.3], opts={"allow_function_pointers": fp}, extralibs=[so],
                         extralibdirs=[os.path.dirname(os.path.abspath(so))])
        flag, outs, S, needed = self.call(lib, g, [0.3])
        self.assertEqual(flag, 0)
        self.checkarray(outs[0], np.sin(0.3) + np.sin(0.6))
        check(S.to_native(), lib_fp if fp else lib_stats, name)
      # VM, also after deserialization: only libraries with function pointers
      for h in [g, ca.Function.deserialize(g.serialize())]:
        V = ca.StatsRecorder({"size": 4096})
        self.checkarray(h(V, 0.3), np.sin(0.3) + np.sin(0.6))
        self.assertFalse(V.truncated())
        check(V.to_native(), lib_fp, name)
      # Found through <name>_config_root, also under a compile-time prefix
      rn = "frenamed_" + name
      fr = ca.Function(rn, [x], [ca.sin(x)])
      cg = ca.CodeGenerator(rn, lib_opts)
      cg.add(fr)
      cg.generate()
      Fr, rso = self.compile_external(rn, rn + ".c",
                                      definitions=["CASADI_CODEGEN_PREFIX", "CODEGEN_PREFIX=zz_"])
      config_root = getattr(cdll(rso, "zz"), rn + "_config_root")
      config_root.restype = ctypes.c_char_p
      self.assertEqual(config_root(), b"zz_")
      g = ca.Function("g", [y], [Fr(y) + Fr(2 * y)])
      for h in [g, ca.Function.deserialize(g.serialize())]:
        V = ca.StatsRecorder({"size": 4096})
        self.checkarray(h(V, 0.3), np.sin(0.3) + np.sin(0.6))
        self.assertFalse(V.truncated())
        check(V.to_native(), lib_fp, "zz_", rn)

  def test_codegen_function_pointers(self):
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    for fp in [True, False]:
      cg = ca.CodeGenerator("fp", {"stats": True, "with_header": True,
                                   "allow_function_pointers": fp})
      cg.add(f)
      with tempfile.TemporaryDirectory() as d:
        cg.generate(d + os.sep)
        with open(os.path.join(d, "fp.c")) as c: code = c.read()
        with open(os.path.join(d, "fp.h")) as h: header = h.read()
      self.assertIn("int CASADI_PREFIX(stats_function_pointers)(void) { return %d; }"
                    % (1 if fp else 0), code)
      self.assertIn("int fp_stats_function_pointers(void);", header)
      self.assertIn("const char* f_config_root(void) { return CASADI_PREFIX_STRING; }", code)
      self.assertNotIn("f_stats_function_pointers", code + header)
      self.assertEqual("(*reserve)" in code, fp)
      self.assertEqual("return s->reserve(s, len, pos);" in code, fp)
      self.assertEqual("if (s->reserve || s->commit) return 0;" in code, not fp)
      self.assertNotIn("#if", code[code.index("casadi_stats_reserve(struct"):][:600])
      self.assertIn("  unsigned char* p;\n  casadi_int cap;\n  casadi_int needed;\n", header)
      self.assertEqual("  void* reserve;\n  void* commit;\n  void* data;\n" in header, not fp)
      self.assertEqual("(*commit)" in header, fp)
      # No mutex without thread_safe
      self.assertIn("*pos = s->needed; s->needed += len;", code)
      self.assertNotIn("mutex", code + header)
      self.assertIn("CASADI_SYMBOL_EXPORT unsigned char* casadi_stats_reserve(struct casadi_stats_sink* s",
                    code)
      self.assertIn("unsigned char* fp_stats_reserve(struct casadi_stats_sink* s, casadi_int len, "
                    "casadi_int* pos);", header)
      self.assertIn("struct casadi_stats_sink fp_stats_make_sink(unsigned char* p, casadi_int cap);",
                    header)
      self.assertIn("CASADI_SYMBOL_EXPORT struct casadi_stats_sink casadi_stats_make_sink(unsigned char* p", code)
      self.assertIn("void fp_stats_clear(struct casadi_stats_sink* s);", header)
      self.assertIn("casadi_int fp_stats_nbytes(const struct casadi_stats_sink* s);", header)
      self.assertIn("int fp_stats_truncated(const struct casadi_stats_sink* s);", header)
      self.assertIn("const unsigned char* fp_stats_data(const struct casadi_stats_sink* s, "
                    "casadi_int* n);", header)
      self.assertNotIn("#define casadi_stats_sink ", code)
      for kind in ["int", "real", "text"]:
        self.assertIn("CASADI_SYMBOL_EXPORT int casadi_get_stat_%s(const struct casadi_stats_sink* s"
                      % kind, code)
        self.assertIn("int fp_get_stat_%s(const struct casadi_stats_sink* s" % kind, header)
      self.assertNotIn("CASADI_MUTEX", code)
      self.assertNotIn("CASADI_ATOMIC", code)
      self.assertNotIn("pthread", code)
      self.assertNotIn("thread_type", code)
      self.assertIn("struct casadi_stats_sink {", header)
      self.assertIn("fp_get_stat_real(", header)
      self.assertNotIn("CASADI_MUTEX_TYPE", header)
      self.assertNotIn("pthread", header)

  @requiresPlugin(ca.Importer, "shell")
  @slow()
  def test_jit(self):
    # Call tree continues into jit-compiled code, also after deserialization
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y")
    G = ca.Function("g", [y], [f(y) + f(2 * y)])
    V = ca.StatsRecorder()
    G(V, 0.3)
    full = shape(V.to_native()[0])
    for codegen_options, inside in [({}, False), ({"stats": True}, True),
                                    ({"stats": True, "allow_function_pointers": False}, False)]:
      g = ca.Function("g", [y], [f(y) + f(2 * y)],
                      {"jit": True, "compiler": "shell", "codegen_options": codegen_options})
      for h in [g, ca.Function.deserialize(g.serialize())]:
        V = ca.StatsRecorder()
        self.checkarray(h(V, 0.3), G(0.3))
        root = V.to_native()[0]
        if inside:
          [c] = root["children"]
          self.assertTrue(c["id"].startswith("jit"))
          self.assertEqual(shape(c), full)
        else:
          self.assertEqual(root["children"], [])

  @requires_integrator("cvodes")
  @requiresPlugin(ca.Importer, "shell")
  @slow()
  def test_jit_oracle(self):
    # Call tree continues into jit-compiled oracle functions
    x = ca.MX.sym("x")
    dae = {"x": x, "ode": -x}
    for codegen_options, inside in [({}, False), ({"stats": True}, True)]:
      I = ca.integrator("I", "cvodes", dae, 0, 1,
                        {"jit": True, "compiler": "shell", "codegen_options": codegen_options})
      V = ca.StatsRecorder()
      self.checkarray(I(V, x0=1)["xf"], np.exp(-1), digits=5)
      calls = [c for c in find(V.to_native(), "daeF") if c["id"].startswith("#")]
      self.assertTrue(calls)
      for c in calls:
        self.assertEqual([d["name"] for d in c["children"]], ["daeF"] if inside else [])

  def test_codegen_parallel_map_needs_thread_safe(self):
    # Without thread_safe, the counter is not atomic
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y", 1, 4)
    for par in ["openmp", "thread"]:
      g = ca.Function("g", [y], [f.map(4, par)(y)])
      cg = ca.CodeGenerator("g", {"stats": True})
      with self.assertInException("requires option thread_safe"):
        cg.add(g)
      cg = ca.CodeGenerator("g", {"stats": True, "thread_safe": True})
      cg.add(g)
      with tempfile.TemporaryDirectory() as d:
        cg.generate(d + os.sep)
        with open(os.path.join(d, "g.c")) as c: code = c.read()
      self.assertIn("#if defined(CASADI_ATOMIC_FETCH_ADD)\n  *pos = CASADI_ATOMIC_FETCH_ADD(&s->needed, len);\n"
                    "#elif defined(CASADI_THREAD_TYPE)\n  CASADI_MUTEX_LOCK(&casadi_stats_mutex);\n", code)

  @slow()
  def test_codegen_stats_off(self):
    # Without the option, generated code does not change
    x = ca.SX.sym("x")
    f = ca.Function("f", [x], [ca.sin(x)])
    y = ca.MX.sym("y", 1, 4)
    g = ca.Function("g", [y], [f.map(4, "thread")(y)])
    a = ca.CodeGenerator("a", {}); a.add(g)
    b = ca.CodeGenerator("a", {"stats": False}); b.add(g)
    self.assertEqual(a.dump(), b.dump())
    self.assertFalse("stats" in a.dump())
    # No file-scope reference counting either (#4324)
    self.assertNotIn("casadi_file_", a.dump())
    # With the option, the plain entry point still works (no sink)
    flags, defs = parallel_options("thread")
    self.check_codegen(g, inputs=[np.arange(4.)], opts={"stats": True, "thread_safe": True},
                       std="c99", extra_options=flags, definitions=defs)


if __name__ == '__main__':
    unittest.main()
