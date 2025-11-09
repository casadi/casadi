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
#
# Tests for the GraphBuilder facade: numeric black-box evaluation through ONNX Runtime,
# and symbolic ONNX import/export. Consolidates the former onnx_runtime.py + onnx_symbolic.py.
import casadi as ca
import numpy
import unittest
import os
import tempfile
from helpers import *

# onnx python package: builds the numeric test models
try:
  import onnx  # type: ignore
  from onnx import helper, TensorProto, numpy_helper  # type: ignore
  have_onnx = True
except ImportError:
  have_onnx = False

# onnxruntime python package: cross-checks the symbolic export
try:
  import onnxruntime as ort  # type: ignore
  has_onnxruntime = True
except ImportError:
  has_onnxruntime = False

# the onnxruntime numeric backend (built alongside the GraphModel onnx backend)
try:
  have_backend = ca.has_onnx("ort")
except Exception:
  have_backend = False

# has_onnx("ort") is True even with the bundled mockup stub (which yields a NULL OrtApi).
# The real ONNX Runtime is supplied on the loader path like every other mockup (commercial_solvers;
# the plugin's RUNPATH=$ORIGIN lets LD_LIBRARY_PATH win over the shipped stub). Probe an actual
# create to tell real from mockup, and gate the numeric suite on it.
have_real_ort = False
if have_onnx and have_backend:
  try:
    import tempfile as _tf
    _fd, _probe = _tf.mkstemp(suffix=".onnx"); os.close(_fd)
    _g = helper.make_graph(
        [helper.make_node("Add", ["x", "a"], ["y"])], "probe",
        [helper.make_tensor_value_info("x", TensorProto.FLOAT, [1])],
        [helper.make_tensor_value_info("y", TensorProto.FLOAT, [1])],
        [numpy_helper.from_array(numpy.full((1,), 1.0, numpy.float32), "a")])
    _m = helper.make_model(_g, opset_imports=[helper.make_opsetid("", 13)]); _m.ir_version = 8
    onnx.save(_m, _probe)
    ca.GraphBuilder(_probe).create("ort_probe")   # raises with the mockup (NULL OrtApi)
    have_real_ort = True
  except Exception:
    have_real_ort = False
  finally:
    try:
      os.remove(_probe)
    except Exception:
      pass

# the GraphModel onnx (symbolic import/export) backend; absent on WITH_ONNX=OFF builds
try:
  have_graph_onnx = ca.has_graphmodel("onnx")
except Exception:
  have_graph_onnx = False

# ONNX element type per numpy dtype, for building test models
ELEM = {numpy.float32: TensorProto.FLOAT, numpy.float64: TensorProto.DOUBLE} if have_onnx else {}

# ONNX Runtime has no float64 kernel for these ops, so we validate them via a float32
# export (casadi_real="float") at reduced precision instead of the default float64 path.
ort_float32 = {
  "unary_cos", "unary_tan", "unary_asin", "unary_acos", "unary_atan",
  "unary_sinh", "unary_cosh", "unary_asinh", "unary_acosh", "unary_atanh", "unary_erf",
  "atan2", "det",
}

# ONNX Runtime can't validate these at all (0-input/shape semantics); roundtrip still runs.
ort_skip = {"empty_function"}

# Cases validated by round-trip only (ORT execution not checked).
ort_roundtrip_only = set()

# Unary operations: (name, expression, list of test inputs)
unary_ops = [
  ("sin", lambda x: ca.sin(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-0.5)]),
  ("cos", lambda x: ca.cos(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-0.5)]),
  ("tan", lambda x: ca.tan(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0)]),
  ("exp", lambda x: ca.exp(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-1.0)]),
  ("log", lambda x: ca.log(x), [ca.DM(0.1), ca.DM(1.0), ca.DM(2.0), ca.DM(10.0)]),
  ("sqrt", lambda x: ca.sqrt(x), [ca.DM(0.0), ca.DM(1.0), ca.DM(4.0), ca.DM(9.0)]),
  ("asin", lambda x: ca.asin(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-0.5)]),
  ("acos", lambda x: ca.acos(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0)]),
  ("atan", lambda x: ca.atan(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-1.0)]),
  ("sinh", lambda x: ca.sinh(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-0.5)]),
  ("cosh", lambda x: ca.cosh(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0)]),
  ("tanh", lambda x: ca.tanh(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-1.0)]),
  ("asinh", lambda x: ca.asinh(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-0.5)]),
  ("acosh", lambda x: ca.acosh(x), [ca.DM(1.0), ca.DM(1.5), ca.DM(2.0)]),
  ("atanh", lambda x: ca.atanh(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(-0.5)]),
  ("ceil", lambda x: ca.ceil(x), [ca.DM(1.2), ca.DM(1.5), ca.DM(1.9), ca.DM(-1.2)]),
  ("floor", lambda x: ca.floor(x), [ca.DM(1.2), ca.DM(1.5), ca.DM(1.9), ca.DM(-1.2)]),
  ("fabs", lambda x: ca.fabs(x), [ca.DM(1.0), ca.DM(-1.0), ca.DM(0.0), ca.DM(-5.5)]),
  ("sign", lambda x: ca.sign(x), [ca.DM(1.0), ca.DM(-1.0), ca.DM(0.0), ca.DM(5.5)]),
  ("neg", lambda x: -x, [ca.DM(1.0), ca.DM(-1.0), ca.DM(0.0), ca.DM(5.5)]),
  ("erf", lambda x: ca.erf(x), [ca.DM(0.0), ca.DM(0.5), ca.DM(1.0), ca.DM(-0.5)]),
]

# Binary operations: (name, expression, list of (x, y) test inputs)
binary_ops = [
  ("add", lambda x, y: x + y, [(ca.DM(2.0), ca.DM(3.0)), (ca.DM(-1.0), ca.DM(4.0)), (ca.DM(0.0), ca.DM(0.0))]),
  ("sub", lambda x, y: x - y, [(ca.DM(5.5), ca.DM(2.2)), (ca.DM(10.0), ca.DM(3.0)), (ca.DM(0.0), ca.DM(0.0))]),
  ("mul", lambda x, y: x * y, [(ca.DM(2.0), ca.DM(3.0)), (ca.DM(-2.0), ca.DM(4.0)), (ca.DM(0.0), ca.DM(5.0))]),
  ("div", lambda x, y: x / y, [(ca.DM(10.0), ca.DM(2.0)), (ca.DM(7.0), ca.DM(3.0)), (ca.DM(1.0), ca.DM(1.0))]),
  ("pow", lambda x, y: x ** y, [(ca.DM(2.0), ca.DM(3.0)), (ca.DM(4.0), ca.DM(0.5)), (ca.DM(10.0), ca.DM(2.0))]),
]


# ============================ Numeric black-box suite ============================
@unittest.skipUnless(have_onnx and have_backend and have_real_ort,
                     "needs the onnx python package and a real (preloaded) ONNX Runtime")
class GraphBuilderNumericTests(casadiTestCase):

  def save_model(self, graph):
    # Serialize a GraphProto to a temporary .onnx file and return its path
    model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
    model.ir_version = 8
    onnx.checker.check_model(model)
    fd, path = tempfile.mkstemp(suffix=".onnx")
    os.close(fd)
    onnx.save(model, path)
    self._tmp.append(path)
    return path

  def setUp(self):
    self._tmp = []

  def tearDown(self):
    for p in self._tmp:
      if os.path.exists(p):
        os.remove(p)

  def affine_model(self, dtype, shape, a, b):
    # y = a*x + b, elementwise, over a tensor of the given shape and dtype
    et = ELEM[dtype]
    av = numpy_helper.from_array(numpy.full(shape, a, dtype), "a")
    bv = numpy_helper.from_array(numpy.full(shape, b, dtype), "b")
    x = helper.make_tensor_value_info("x", et, list(shape))
    y = helper.make_tensor_value_info("y", et, list(shape))
    g = helper.make_graph(
        [helper.make_node("Mul", ["x", "a"], ["t"]),
         helper.make_node("Add", ["t", "b"], ["y"])],
        "affine", [x], [y], [av, bv])
    return self.save_model(g)

  def test_vector(self):
    self.message("black-box eval of a float32 vector model")
    path = self.affine_model(numpy.float32, (3,), 2.0, 1.0)
    f = ca.GraphBuilder(path).create("f")
    self.assertEqual(f.n_in(), 1)
    self.assertEqual(f.n_out(), 1)
    self.assertEqual(f.name_in(), ["x"])
    self.assertEqual(f.name_out(), ["y"])
    self.checkarray(f.size_in(0), (3, 1), "vector maps to column")
    x = ca.DM([0.5, 1.0, -2.0])
    self.checkarray(f(x), 2.0 * x + 1.0, "vector", digits=5)

  def test_double(self):
    self.message("black-box eval of a float64 (double) model")
    path = self.affine_model(numpy.float64, (4,), -3.0, 0.5)
    f = ca.GraphBuilder(path).create("f")
    x = ca.DM([1.0, 2.0, 3.0, 4.0])
    self.checkarray(f(x), -3.0 * x + 0.5, "double", digits=12)

  def test_serialize(self):
    self.message("an OnnxFunction round-trips through serialization")
    path = self.affine_model(numpy.float64, (4,), -3.0, 0.5)
    f = ca.GraphBuilder(path).create("f")
    self.check_serialize(f, [ca.DM([1.0, 2.0, 3.0, 4.0])])

  def test_matrix_layout(self):
    self.message("rank-2 tensor: ONNX row-major vs CasADi column-major reconciliation")
    path = self.affine_model(numpy.float32, (2, 3), 1.0, 0.0)  # identity-ish: y = x
    f = ca.GraphBuilder(path).create("f")
    self.checkarray(f.size_in(0), (2, 3), "2-D stays 2-D")
    x = ca.DM([[1, 2, 3], [4, 5, 6]])
    self.checkarray(f(x), x, "matrix is not transposed", digits=5)

  def test_highrank_flatten(self):
    self.message("rank-3 tensor flattens to a column vector")
    path = self.affine_model(numpy.float32, (2, 2, 2), 2.0, 0.0)
    f = ca.GraphBuilder(path).create("f")
    self.checkarray(f.size_in(0), (8, 1), "rank>2 flattens")
    x = ca.DM(range(8))
    self.checkarray(f(x), 2.0 * x, "flattened", digits=5)

  def test_highrank_order(self):
    self.message("rank-4 tensor: create() succeeds; row-major flatten order is preserved")
    shape = (2, 3, 2, 2)  # 24 elements, rank 4
    n = int(numpy.prod(shape))
    # Distinct value per element so a wrong flatten order would mismatch (unlike elementwise *2)
    c = numpy.arange(n, dtype=numpy.float32).reshape(shape)
    et = ELEM[numpy.float32]
    x = helper.make_tensor_value_info("x", et, list(shape))
    y = helper.make_tensor_value_info("y", et, list(shape))
    g = helper.make_graph(
        [helper.make_node("Add", ["x", "c"], ["y"])],
        "highrank", [x], [y], [numpy_helper.from_array(c, "c")])
    f = ca.GraphBuilder(self.save_model(g)).create("f")
    self.checkarray(f.size_in(0), (n, 1), "rank>2 input flattens to a column")
    self.checkarray(f.size_out(0), (n, 1), "rank>2 output flattens to a column")
    xv = ca.DM(range(100, 100 + n))
    # ONNX is row-major (C-order); the flattened column must follow that ordering
    expected = xv + ca.DM(c.flatten(order="C").tolist())
    self.checkarray(f(xv), expected, "rank-4 add, row-major order", digits=4)

  def test_multi_io(self):
    self.message("two inputs, two outputs; subset selection by name")
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, [2])
    y = helper.make_tensor_value_info("y", TensorProto.FLOAT, [2])
    s = helper.make_tensor_value_info("s", TensorProto.FLOAT, [2])
    d = helper.make_tensor_value_info("d", TensorProto.FLOAT, [2])
    g = helper.make_graph(
        [helper.make_node("Add", ["x", "y"], ["s"]),
         helper.make_node("Sub", ["x", "y"], ["d"])],
        "two", [x, y], [s, d])
    path = self.save_model(g)

    f = ca.GraphBuilder(path).create("f")
    self.assertEqual(f.name_out(), ["s", "d"])
    xv, yv = ca.DM([1, 2]), ca.DM([3, 5])
    r = f(xv, yv)
    self.checkarray(r[0], xv + yv, "sum", digits=5)
    self.checkarray(r[1], xv - yv, "diff", digits=5)

    # Select a single output and reorder inputs
    f2 = ca.GraphBuilder(path).create("f2", ["y", "x"], ["d"])
    self.assertEqual(f2.name_in(), ["y", "x"])
    self.assertEqual(f2.name_out(), ["d"])
    self.checkarray(f2(yv, xv), xv - yv, "selected diff", digits=5)

  def two_pathway_model(self):
    # Two fully independent pathways: ya = a*2 (needs only a), yb = b+10 (needs only b)
    a = helper.make_tensor_value_info("a", TensorProto.FLOAT, [2])
    b = helper.make_tensor_value_info("b", TensorProto.FLOAT, [2])
    ya = helper.make_tensor_value_info("ya", TensorProto.FLOAT, [2])
    yb = helper.make_tensor_value_info("yb", TensorProto.FLOAT, [2])
    two = numpy_helper.from_array(numpy.array([2, 2], numpy.float32), "two")
    ten = numpy_helper.from_array(numpy.array([10, 10], numpy.float32), "ten")
    g = helper.make_graph(
        [helper.make_node("Mul", ["a", "two"], ["ya"]),
         helper.make_node("Add", ["b", "ten"], ["yb"])],
        "two_pathway", [a, b], [ya, yb], [two, ten])
    return self.save_model(g)

  def test_integer_dtype(self):
    self.message("non-float numeric dtypes are converted to/from double")
    one = numpy_helper.from_array(numpy.array([1, 1, 1], numpy.int64), "one")
    x = helper.make_tensor_value_info("x", TensorProto.INT64, [3])
    y = helper.make_tensor_value_info("y", TensorProto.INT64, [3])
    g = helper.make_graph([helper.make_node("Add", ["x", "one"], ["y"])],
                          "int64", [x], [y], [one])
    path = self.save_model(g)
    f = ca.GraphBuilder(path).create("f")
    self.checkarray(f(ca.DM([10, 20, 30])), ca.DM([11, 21, 31]), "int64", digits=12)

  def fwd_model(self, n):
    # One model holding the primal y=2x and its forward sensitivity fwd_y=2*fwd_x, named per
    # CasADi's forward convention. The seed count is the symbolic dimension "nfwd".
    two1 = numpy_helper.from_array(numpy.full((n,), 2, numpy.float32), "two1")
    two2 = numpy_helper.from_array(numpy.full((n, 1), 2, numpy.float32), "two2")
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, [n])
    out_y = helper.make_tensor_value_info("out_y", TensorProto.FLOAT, [n])  # nondiff output (unused)
    fwd_x = helper.make_tensor_value_info("fwd_x", TensorProto.FLOAT, [n, "nfwd"])
    y = helper.make_tensor_value_info("y", TensorProto.FLOAT, [n])
    fwd_y = helper.make_tensor_value_info("fwd_y", TensorProto.FLOAT, [n, "nfwd"])
    g = helper.make_graph(
        [helper.make_node("Mul", ["x", "two1"], ["y"]),
         helper.make_node("Mul", ["fwd_x", "two2"], ["fwd_y"])],
        "withfwd", [x, out_y, fwd_x], [y, fwd_y], [two1, two2])
    return self.save_model(g)

  def test_get_forward(self):
    self.message("a model authored with CasADi fwd naming serves its own forward derivatives")
    path = self.fwd_model(3)
    f = ca.GraphBuilder(path).create("f", ["x"], ["y"])
    self.checkarray(f(ca.DM([1, 2, 3])), ca.DM([2, 4, 6]), "primal", digits=5)

    # Forward function signature: [x, out_y, fwd_x] -> [fwd_y]
    fwd = f.forward(2)
    self.assertEqual(list(fwd.name_in()), ["x", "out_y", "fwd_x"])
    self.assertEqual(list(fwd.name_out()), ["fwd_y"])
    seeds = ca.DM([[1, 4], [2, 5], [3, 6]])
    fy = fwd(ca.DM([1, 2, 3]), ca.DM([2, 4, 6]), seeds)
    self.checkarray(fy, 2 * seeds, "forward sensitivity", digits=5)

    # Jacobian via forward mode is exact (would be noisy if it fell back to finite differences)
    X = ca.MX.sym("X", 3)
    J = ca.Function("J", [X], [ca.jacobian(f(X), X)])
    self.checkarray(J(ca.DM([1, 2, 3])), 2 * ca.DM.eye(3), "forward-mode jacobian", digits=10)

  def rev_model(self, n):
    # Primal y=2x and adjoint adj_x=2*adj_y in one model; seed count is symbolic "nadj"
    two1 = numpy_helper.from_array(numpy.full((n,), 2, numpy.float32), "two1")
    two2 = numpy_helper.from_array(numpy.full((n, 1), 2, numpy.float32), "two2")
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, [n])
    out_y = helper.make_tensor_value_info("out_y", TensorProto.FLOAT, [n])
    adj_y = helper.make_tensor_value_info("adj_y", TensorProto.FLOAT, [n, "nadj"])
    y = helper.make_tensor_value_info("y", TensorProto.FLOAT, [n])
    adj_x = helper.make_tensor_value_info("adj_x", TensorProto.FLOAT, [n, "nadj"])
    g = helper.make_graph(
        [helper.make_node("Mul", ["x", "two1"], ["y"]),
         helper.make_node("Mul", ["adj_y", "two2"], ["adj_x"])],
        "withadj", [x, out_y, adj_y], [y, adj_x], [two1, two2])
    return self.save_model(g)

  def jac_model(self, n):
    # Primal y=2x and its constant jacobian jac_y_x = 2*I in one model
    two1 = numpy_helper.from_array(numpy.full((n,), 2, numpy.float32), "two1")
    eye2 = numpy_helper.from_array((2 * numpy.eye(n)).astype(numpy.float32), "eye2")
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, [n])
    out_y = helper.make_tensor_value_info("out_y", TensorProto.FLOAT, [n])
    y = helper.make_tensor_value_info("y", TensorProto.FLOAT, [n])
    jac = helper.make_tensor_value_info("jac_y_x", TensorProto.FLOAT, [n, n])
    g = helper.make_graph(
        [helper.make_node("Mul", ["x", "two1"], ["y"]),
         helper.make_node("Identity", ["eye2"], ["jac_y_x"])],
        "withjac", [x, out_y], [y, jac], [two1, eye2])
    return self.save_model(g)

  def test_get_reverse(self):
    self.message("a model authored with CasADi adj naming serves its own reverse derivatives")
    f = ca.GraphBuilder(self.rev_model(3)).create("f", ["x"], ["y"])
    rev = f.reverse(2)
    self.assertEqual(list(rev.name_in()), ["x", "out_y", "adj_y"])
    self.assertEqual(list(rev.name_out()), ["adj_x"])
    seeds = ca.DM([[1, 4], [2, 5], [3, 6]])
    self.checkarray(rev(ca.DM([1, 2, 3]), ca.DM([2, 4, 6]), seeds), 2 * seeds,
                    "adjoint sensitivity", digits=5)

  def test_get_jacobian(self):
    self.message("a model authored with CasADi jac naming serves its own jacobian")
    f = ca.GraphBuilder(self.jac_model(3)).create("f", ["x"], ["y"])
    J = f.jacobian()
    self.assertEqual(list(J.name_in()), ["x", "out_y"])
    self.assertEqual(list(J.name_out()), ["jac_y_x"])
    self.checkarray(J(ca.DM([1, 2, 3]), ca.DM([2, 4, 6])), 2 * ca.DM.eye(3), "jacobian", digits=8)
    # ca.jacobian routes through get_jacobian -> exact
    X = ca.MX.sym("X", 3)
    Jf = ca.Function("Jf", [X], [ca.jacobian(f(X), X)])
    self.checkarray(Jf(ca.DM([1, 2, 3])), 2 * ca.DM.eye(3), "jacobian via graph", digits=10)

  def isdiff_fwd_model(self, n):
    # y = 2x + p, with forward fwd_y = 2*fwd_x. p is meant to be non-differentiable, so the
    # model provides NO fwd_p seed and NO out_y (the wrapper must tolerate their absence).
    two1 = numpy_helper.from_array(numpy.full((n,), 2, numpy.float32), "two1")
    two2 = numpy_helper.from_array(numpy.full((n, 1), 2, numpy.float32), "two2")
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, [n])
    p = helper.make_tensor_value_info("p", TensorProto.FLOAT, [n])
    fwd_x = helper.make_tensor_value_info("fwd_x", TensorProto.FLOAT, [n, "nfwd"])
    y = helper.make_tensor_value_info("y", TensorProto.FLOAT, [n])
    fwd_y = helper.make_tensor_value_info("fwd_y", TensorProto.FLOAT, [n, "nfwd"])
    g = helper.make_graph(
        [helper.make_node("Mul", ["x", "two1"], ["t"]),
         helper.make_node("Add", ["t", "p"], ["y"]),
         helper.make_node("Mul", ["fwd_x", "two2"], ["fwd_y"])],
        "isdiff", [x, p, fwd_x], [y, fwd_y], [two1, two2])
    return self.save_model(g)

  def test_is_diff_forward(self):
    self.message("is_diff_in=false inputs need no fwd_ seed; absent out_/fwd_ are tolerated")
    path = self.isdiff_fwd_model(3)
    f = ca.GraphBuilder(path).create("f", ["x", "p"], ["y"],
        {"is_diff_in": [True, False]})
    self.checkarray(f(ca.DM([1, 2, 3]), ca.DM([10, 20, 30])), ca.DM([12, 24, 36]),
                    "primal", digits=5)

    # Succeeds without enable_fd, so the analytic forward (not finite differences) is used
    fwd = f.forward(2)
    self.assertEqual(list(fwd.name_in()), ["x", "p", "out_y", "fwd_x", "fwd_p"])
    self.assertEqual(list(fwd.name_out()), ["fwd_y"])
    seeds = ca.DM([[1, 4], [2, 5], [3, 6]])
    fy = fwd(ca.DM([1, 2, 3]), ca.DM([10, 20, 30]), ca.DM([12, 24, 36]),
             seeds, ca.DM.zeros(3, 2))
    self.checkarray(fy, 2 * seeds, "forward ignores the non-diff seed", digits=5)

    # Jacobian w.r.t the differentiable input, via forward mode, is exact
    X = ca.MX.sym("x", 3)
    P = ca.MX.sym("p", 3)
    J = ca.Function("J", [X, P], [ca.jacobian(f(X, P), X)])
    self.checkarray(J(ca.DM([1, 2, 3]), ca.DM([10, 20, 30])), 2 * ca.DM.eye(3),
                    "jacobian wrt x", digits=10)

  def test_partial_inference(self):
    self.message("select one independent pathway; ORT runs only that subgraph")
    path = self.two_pathway_model()

    # Wire through pathway A only: needs only input a, produces only ya
    fa = ca.GraphBuilder(path).create("fa", ["a"], ["ya"])
    self.assertEqual(list(fa.name_in()), ["a"])
    self.assertEqual(list(fa.name_out()), ["ya"])
    self.checkarray(fa(ca.DM([1, 2])), ca.DM([2, 4]), "pathway A", digits=5)

    # Wire through pathway B only: needs only input b, produces only yb
    fb = ca.GraphBuilder(path).create("fb", ["b"], ["yb"])
    self.checkarray(fb(ca.DM([1, 2])), ca.DM([11, 12]), "pathway B", digits=5)

    # The full model still works with both pathways
    f = ca.GraphBuilder(path).create("f")
    r = f(ca.DM([1, 2]), ca.DM([3, 4]))
    self.checkarray(r[0], ca.DM([2, 4]), "full ya", digits=5)
    self.checkarray(r[1], ca.DM([13, 14]), "full yb", digits=5)

    # Asking for yb without wiring its input b: b defaults to NaN, poisoning the output
    fbad = ca.GraphBuilder(path).create("fbad", ["a"], ["yb"])
    self.assertTrue(numpy.all(numpy.isnan(numpy.array(fbad(ca.DM([1, 2]))))))

  def test_finite_differences(self):
    self.message("derivatives of a black box fall back to finite differences")
    path = self.affine_model(numpy.float64, (3,), 2.0, 1.0)
    f = ca.GraphBuilder(path).create("f", {"enable_fd": True})
    x = ca.MX.sym("x", 3)
    J = ca.Function("J", [x], [ca.jacobian(f(x), x)])
    self.checkarray(J(ca.DM([0.5, 1.0, -2.0])), 2.0 * ca.DM.eye(3), "fd jacobian", digits=4)

  def dyn_model(self):
    # y = x*2, input "x" with a symbolic batch dimension: [batch, 3]
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, ["batch", 3])
    y = helper.make_tensor_value_info("y", TensorProto.FLOAT, ["batch", 3])
    two = numpy_helper.from_array(numpy.array([[2, 2, 2]], numpy.float32), "two")
    g = helper.make_graph([helper.make_node("Mul", ["x", "two"], ["y"])],
                          "dyn", [x], [y], [two])
    return self.save_model(g)

  def test_builder_introspect(self):
    self.message("GraphBuilder model introspection")
    path = self.affine_model(numpy.float32, (2, 3), 1.0, 0.0)
    b = ca.GraphBuilder(path)
    self.assertEqual(b.n_in(), 1)
    self.assertEqual(b.n_out(), 1)
    self.assertEqual(list(b.name_in()), ["x"])
    self.assertEqual(list(b.name_out()), ["y"])
    self.assertEqual(list(b.input_shape("x")), [2, 3])
    self.assertEqual(list(b.dynamic_params()), [])

  def test_builder_dynamic_dim(self):
    self.message("GraphBuilder dynamic-dimension binding")
    path = self.dyn_model()
    b = ca.GraphBuilder(path)
    self.assertEqual(list(b.dynamic_params()), ["batch"])
    self.assertEqual(list(b.input_shape("x")), [-1, 3])
    b.bind_dim("batch", 2)
    f = b.create("net", b.name_in(), b.name_out())
    self.checkarray(f.size_in(0), (2, 3), "bound batch")
    x = ca.DM([[1, 2, 3], [4, 5, 6]])
    self.checkarray(f(x), 2.0 * x, "dynamic eval", digits=5)

  def test_codegen_smoke(self):
    self.message("C code generation embeds the model and calls the runtime")
    path = self.affine_model(numpy.float32, (3,), 2.0, 1.0)
    f = ca.GraphBuilder(path).create("f")
    d = tempfile.mkdtemp()
    cwd = os.getcwd()
    try:
      os.chdir(d)
      f.generate("f.c")
      with open(os.path.join(d, "f.c")) as fh:
        code = fh.read()
    finally:
      os.chdir(cwd)
    self.assertTrue("casadi_onnxruntime_solve" in code)
    self.assertTrue("casadi_onnxruntime_prob" in code)
    self.assertTrue("ort_runtime.h" in code)

  def test_codegen_compile(self):
    # Full roundtrip: generate C, compile+link against ONNX Runtime, load and evaluate.
    # Gated on ONNXRUNTIME_DIR (path to an ONNX Runtime C/C++ distribution).
    ort_dir = os.environ.get("ONNXRUNTIME_DIR")
    if not ort_dir:
      self.skipTest("set ONNXRUNTIME_DIR to a C/C++ distribution to run the codegen compile test")
    import subprocess
    self.message("generated C compiles against ONNX Runtime and evaluates correctly")
    path = self.affine_model(numpy.float64, (3,), 2.0, 1.0)
    f = ca.GraphBuilder(path).create("f")
    d = tempfile.mkdtemp()
    cwd = os.getcwd()
    try:
      os.chdir(d)
      f.generate("f.c")
      inc = os.path.join(os.path.dirname(ca.__file__), "include",
                         "casadi", "interfaces", "ort")
      r = subprocess.run(["gcc", "-shared", "-fPIC", "f.c",
                          "-I" + inc, "-I" + os.path.join(ort_dir, "include"),
                          "-L" + os.path.join(ort_dir, "lib"), "-lonnxruntime", "-o", "f.so"],
                         capture_output=True, text=True)
      self.assertEqual(r.returncode, 0, r.stderr)
      g = ca.external("f", os.path.join(d, "f.so"))
      x = ca.DM([0.5, 1.0, -2.0])
      self.checkarray(g(x), 2.0 * x + 1.0, "codegen eval", digits=10)
    finally:
      os.chdir(cwd)

  def test_builder_set(self):
    self.message("GraphBuilder.set bakes an input value, consumed at create()")
    path = self.two_pathway_model()
    b = ca.GraphBuilder(path)
    b.set("b", [5.0, 5.0])
    f = b.create("f", b.name_in(), b.name_out())
    self.assertEqual(list(f.name_in()), ["a"])  # b is baked, no longer an input
    r = f(ca.DM([1, 2]))
    self.checkarray(r[0], ca.DM([2, 4]), "ya", digits=5)
    self.checkarray(r[1], ca.DM([15, 15]), "yb uses baked b", digits=5)

  def test_builder_node(self):
    self.message("GraphBuilder.node exposes tensor metadata")
    path = self.dyn_model()
    b = ca.GraphBuilder(path)
    n = b.node("x")
    self.assertEqual(n.name, "x")
    self.assertEqual(n.io, "input")
    self.assertEqual(n.dtype, "FLOAT")
    self.assertEqual(list(n.shape), [-1, 3])
    self.assertEqual(b.node("y").io, "output")

  def test_builder_print(self):
    self.message("GraphBuilder and Node printing")
    b = ca.GraphBuilder(self.two_pathway_model())
    self.assertTrue("GraphBuilder" in str(b) and "input(s)" in str(b))
    n = b.node("a")
    self.assertTrue("FLOAT" in str(n))
    self.assertEqual(n.io, "input")

  def test_builder_bind_shape(self):
    self.message("GraphBuilder explicit shape binding")
    path = self.dyn_model()
    b = ca.GraphBuilder(path)
    b.bind_shape("x", [4, 3])
    f = b.create("net", b.name_in(), b.name_out())
    self.checkarray(f.size_in(0), (4, 3), "pinned shape")
    x = ca.DM(numpy.arange(12).reshape(4, 3))
    self.checkarray(f(x), 2.0 * x, "bound-shape eval", digits=5)



# ============================ Symbolic import/export suite =======================
@unittest.skipUnless(have_graph_onnx, "needs the GraphModel onnx (symbolic) backend")
class GraphBuilderSymbolicTests(casadiTestCase):

  def compare(self, ref, got, name):
    # A CasADi Function returns a DM for a single output, a list for several, or a
    # dict (keyed by output name) when called with no positional arguments
    if isinstance(ref, dict):
      for k in ref:
        self.checkarray(ref[k], got[k], name)
    elif isinstance(ref, (list, tuple)):
      for a, b in zip(ref, got):
        self.checkarray(a, b, name)
    else:
      self.checkarray(ref, got, name)

  def roundtrip(self, name, f, inputs):
    # CasADi -> ONNX -> CasADi, then check the imported function reproduces f
    fn = "test_%s.onnx" % name
    try:
      ca.GraphBuilder(f).export_onnx(fn)
      g = ca.GraphBuilder(fn).create("imported_%s" % name, {"symbolic": True})
      for vals in inputs:
        if not isinstance(vals, tuple):
          vals = (vals,)
        self.compare(f(*vals), g(*vals), name)
    finally:
      if os.path.exists(fn):
        os.remove(fn)

  def onnxruntime_check(self, name, f, inputs, float32=False, digits=10):
    # Export and run through ONNX Runtime, then check it reproduces f. float32 uses the
    # casadi_real="float" export so ORT ops lacking a float64 kernel can still be validated.
    fn = "test_%s_ort.onnx" % name
    dtype = numpy.float32 if float32 else numpy.float64
    try:
      ca.GraphBuilder(f).export_onnx(fn, {"casadi_real": "float" if float32 else "double"})
      session = ort.InferenceSession(fn)
      input_names = [i.name for i in session.get_inputs()]
      for vals in inputs:
        if not isinstance(vals, tuple):
          vals = (vals,)
        feed = {input_names[i]: numpy.array(v).astype(dtype) for i, v in enumerate(vals)}
        ref = f(*vals)
        if not isinstance(ref, (list, tuple)):
          ref = [ref]
        for a, b in zip(ref, session.run(None, feed)):
          self.checkarray(a, ca.DM(b), name, digits=digits)
    finally:
      if os.path.exists(fn):
        os.remove(fn)

  def check(self, name, f, inputs):
    # Every case round-trips (float64). ORT validation: float64 where supported, else a
    # float32 export at reduced precision; a few cases ORT can't validate at all.
    self.message(name)
    self.roundtrip(name, f, inputs)
    if not has_onnxruntime or name in ort_skip or name in ort_roundtrip_only:
      return
    if name in ort_float32:
      self.onnxruntime_check(name, f, inputs, float32=True, digits=4)
    else:
      self.onnxruntime_check(name, f, inputs)

  def test_unary(self):
    self.message("ONNX unary operations")
    x = ca.MX.sym("x")
    for name, expr, inputs in unary_ops:
      self.check("unary_" + name, ca.Function("test_" + name, [x], [expr(x)]), inputs)

  def test_binary(self):
    self.message("ONNX binary operations")
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    for name, expr, inputs in binary_ops:
      self.check("binary_" + name, ca.Function("test_" + name, [x, y], [expr(x, y)]), inputs)

  def test_matrix(self):
    self.message("ONNX matmul and transpose")
    A = ca.MX.sym("A", 2, 3)
    B = ca.MX.sym("B", 3, 2)
    self.check("matmul_2x3_3x2", ca.Function("f", [A, B], [ca.mtimes(A, B)]),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM([[1, 2], [3, 4], [5, 6]]))])

    A = ca.MX.sym("A", 3, 1)
    B = ca.MX.sym("B", 1, 3)
    self.check("matmul_3x1_1x3", ca.Function("f", [A, B], [ca.mtimes(A, B)]),
               [(ca.DM([[1], [2], [3]]), ca.DM([[4, 5, 6]]))])

    A = ca.MX.sym("A", 2, 3)
    self.check("transpose_2x3", ca.Function("f", [A], [A.T]), [ca.DM([[1, 2, 3], [4, 5, 6]])])

  def test_shapes(self):
    self.message("ONNX sin over various input shapes")
    for shape, val in [((1, 1), ca.DM([[2.0]])),
                       ((3, 1), ca.DM([[1.0], [2.0], [3.0]])),
                       ((1, 3), ca.DM([[1.0, 2.0, 3.0]])),
                       ((2, 2), ca.DM([[1.0, 2.0], [3.0, 4.0]]))]:
      x = ca.MX.sym("x", shape[0], shape[1])
      self.check("sin_%dx%d" % shape, ca.Function("f", [x], [ca.sin(x)]), [val])

  def test_expressions(self):
    self.message("ONNX compound expressions and multiple outputs")
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    self.check("empty_function", ca.Function("empty", [], [ca.MX(42.0)]), [()])

    self.check("complex_expression",
               ca.Function("f", [x, y], [(ca.sin(x) + ca.sin(y)) * ca.exp(x - y)]),
               [(ca.DM(0.5), ca.DM(1.0)), (ca.DM(1.0), ca.DM(0.5)), (ca.DM(0.0), ca.DM(0.0))])

    self.check("multiple_outputs",
               ca.Function("f", [x, y], [x + y, x - y, x * y]),
               [(ca.DM(3.0), ca.DM(2.0)), (ca.DM(5.0), ca.DM(1.0))])

  def test_function_hierarchy(self):
    self.message("ONNX nested function calls (OP_CALL)")
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")

    inner = ca.Function("inner", [x], [x + x])
    self.check("function_hierarchy_simple",
               ca.Function("outer", [y], [inner(y) + 1]).wrap(), [ca.DM(5.0)])

    double = ca.Function("double", [x], [x + x])
    self.check("function_hierarchy_multiple_calls",
               ca.Function("outer_multi", [y], [double(y) + double(y + 1)]).wrap(), [ca.DM(3.0)])

    func_c = ca.Function("func_c", [x], [ca.sin(x)])
    func_b = ca.Function("func_b", [y], [func_c(y) + 1])
    self.check("function_hierarchy_deep",
               ca.Function("func_a", [z], [func_b(z) * 2]).wrap(), [ca.DM(0.5)])

  def test_vertcat(self):
    self.message("ONNX vertcat as input and as output")
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    self.check("vertcat_input", ca.Function("f", [ca.vertcat(x, y)], [x + y]), [ca.DM([3.0, 2.0])])
    self.check("vertcat_output", ca.Function("f", [x, y], [ca.vertcat(x, y)]), [(ca.DM(3.0), ca.DM(2.0))])

  def test_mx_ops(self):
    self.message("ONNX repmat, indexing, dot, repsum, blockcat")
    x = ca.MX.sym("x", 2, 2)
    self.check("repmat_1x3", ca.Function("f", [x], [ca.repmat(x, 1, 3)]), [ca.DM([[1, 2], [3, 4]])])

    x = ca.MX.sym("x", 3, 1)
    self.check("repmat_1x2", ca.Function("f", [x], [ca.repmat(x, 1, 2)]), [ca.DM([[1], [2], [3]])])

    x = ca.MX.sym("x", 5, 1)
    self.check("index_single_element", ca.Function("f", [x], [x[0]]), [ca.DM([1, 2, 3, 4, 5])])
    self.check("index_multiple_elements", ca.Function("f", [x], [x[[0, 2, 4]]]), [ca.DM([10, 20, 30, 40, 50])])

    x = ca.MX.sym("x", 6, 1)
    self.check("index_slice", ca.Function("f", [x], [x[1:4]]), [ca.DM([1, 2, 3, 4, 5, 6])])

    x = ca.MX.sym("x", 3, 1)
    y = ca.MX.sym("y", 3, 1)
    self.check("dot_product", ca.Function("f", [x, y], [ca.dot(x, y)]), [(ca.DM([1, 2, 3]), ca.DM([4, 5, 6]))])

    x = ca.MX.sym("x", 2, 6)
    self.check("repsum_1x3", ca.Function("f", [x], [ca.repsum(x, 1, 3)]),
               [ca.DM([[1, 2, 10, 20, 100, 200], [3, 4, 30, 40, 300, 400]])])

    a = ca.MX.sym("a", 2, 2)
    b = ca.MX.sym("b", 2, 2)
    c = ca.MX.sym("c", 2, 2)
    d = ca.MX.sym("d", 2, 2)
    self.check("blockcat_2x2", ca.Function("f", [a, b, c, d], [ca.blockcat([[a, b], [c, d]])]),
               [(ca.DM([[1, 2], [3, 4]]), ca.DM([[5, 6], [7, 8]]),
                 ca.DM([[9, 10], [11, 12]]), ca.DM([[13, 14], [15, 16]]))])

  def test_split_concat(self):
    self.message("ONNX horzsplit/vertsplit and horzcat/vertcat")
    x = ca.MX.sym("x", 2, 6)
    parts = ca.horzsplit(x, [0, 2, 4, 6])
    self.check("horzsplit_equal", ca.Function("f", [x], [parts[0] + parts[1] + parts[2]]),
               [ca.DM([[1, 2, 10, 20, 100, 200], [3, 4, 30, 40, 300, 400]])])

    x = ca.MX.sym("x", 3, 5)
    parts = ca.horzsplit(x, [0, 2, 5])
    self.check("horzsplit_unequal", ca.Function("f", [x], [parts[0], parts[1]]),
               [ca.DM([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]])])

    x = ca.MX.sym("x", 6, 2)
    parts = ca.vertsplit(x, [0, 2, 4, 6])
    self.check("vertsplit_equal", ca.Function("f", [x], [parts[0] + parts[1] + parts[2]]),
               [ca.DM([[1, 2], [3, 4], [10, 20], [30, 40], [100, 200], [300, 400]])])

    x = ca.MX.sym("x", 5, 2)
    parts = ca.vertsplit(x, [0, 2, 5])
    self.check("vertsplit_unequal", ca.Function("f", [x], [parts[0], parts[1]]),
               [ca.DM([[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]])])

    a = ca.MX.sym("a", 2, 2)
    b = ca.MX.sym("b", 2, 3)
    self.check("horzcat_2x2_2x3", ca.Function("f", [a, b], [ca.horzcat(a, b)]),
               [(ca.DM([[1, 2], [3, 4]]), ca.DM([[5, 6, 7], [8, 9, 10]]))])

    a = ca.MX.sym("a", 2, 3)
    b = ca.MX.sym("b", 3, 3)
    self.check("vertcat_2x3_3x3", ca.Function("f", [a, b], [ca.vertcat(a, b)]),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM([[7, 8, 9], [10, 11, 12], [13, 14, 15]]))])

  def test_slice_reduce(self):
    self.message("ONNX slicing, reductions, reciprocal and norms")
    x = ca.MX.sym("x", 4, 4)
    self.check("slice_2d", ca.Function("f", [x], [x[1:3, 0:2]]),
               [ca.DM([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])])

    x = ca.MX.sym("x", 4, 3)
    self.check("slice_rows", ca.Function("f", [x], [x[1:3, :]]),
               [ca.DM([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])])

    x = ca.MX.sym("x", 3, 2)
    self.check("reduce_min", ca.Function("f", [x], [ca.mmin(x)]), [ca.DM([[5, 2], [1, 8], [3, 4]])])
    self.check("reduce_max", ca.Function("f", [x], [ca.mmax(x)]), [ca.DM([[5, 2], [1, 8], [3, 4]])])

    x = ca.MX.sym("x", 2, 2)
    self.check("reciprocal", ca.Function("f", [x], [1.0 / x]), [ca.DM([[1, 2], [4, 5]])])

    x = ca.MX.sym("x", 3, 1)
    self.check("norm_1", ca.Function("f", [x], [ca.norm_1(x)]), [ca.DM([1, -2, 3])])
    self.check("norm_2", ca.Function("f", [x], [ca.norm_2(x)]), [ca.DM([3, 4, 0])])

  def test_norms_squares(self):
    self.message("ONNX sq, twice, Frobenius/inf norms")
    x = ca.MX.sym("x", 2, 2)
    self.check("sq", ca.Function("f", [x], [x ** 2]), [ca.DM([[1, 2], [3, 4]])])
    self.check("twice", ca.Function("f", [x], [x + x]), [ca.DM([[1, 2], [3, 4]])])

    x = ca.MX.sym("x", 3, 2)
    self.check("norm_fro", ca.Function("f", [x], [ca.norm_fro(x)]), [ca.DM([[1, 2], [3, 4], [5, 6]])])

    x = ca.MX.sym("x", 3, 1)
    self.check("norm_inf", ca.Function("f", [x], [ca.norm_inf(x)]), [ca.DM([1, -5, 3])])

  def test_comparisons(self):
    self.message("ONNX comparison and logical operations")
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    self.check("ne", ca.Function("f", [x, y], [x != y]), [(ca.DM(1.0), ca.DM(1.0)), (ca.DM(1.0), ca.DM(2.0))])
    self.check("eq", ca.Function("f", [x, y], [x == y]), [(ca.DM(1.0), ca.DM(1.0)), (ca.DM(1.0), ca.DM(2.0))])
    self.check("le", ca.Function("f", [x, y], [x <= y]),
               [(ca.DM(1.0), ca.DM(2.0)), (ca.DM(2.0), ca.DM(2.0)), (ca.DM(3.0), ca.DM(2.0))])
    self.check("lt", ca.Function("f", [x, y], [x < y]),
               [(ca.DM(1.0), ca.DM(2.0)), (ca.DM(2.0), ca.DM(2.0)), (ca.DM(3.0), ca.DM(2.0))])

    self.check("logic_not", ca.Function("f", [x], [ca.logic_not(x)]), [ca.DM(0.0), ca.DM(1.0)])
    self.check("logic_or", ca.Function("f", [x, y], [ca.logic_or(x, y)]),
               [(ca.DM(0.0), ca.DM(0.0)), (ca.DM(0.0), ca.DM(1.0)), (ca.DM(1.0), ca.DM(1.0))])
    self.check("logic_and", ca.Function("f", [x, y], [ca.logic_and(x, y)]),
               [(ca.DM(0.0), ca.DM(0.0)), (ca.DM(0.0), ca.DM(1.0)), (ca.DM(1.0), ca.DM(1.0))])

  def test_binary_extra(self):
    self.message("ONNX fmod, copysign, fmin, fmax")
    x = ca.MX.sym("x", 2, 2)
    y = ca.MX.sym("y", 2, 2)
    self.check("fmod", ca.Function("f", [x, y], [ca.fmod(x, y)]),
               [(ca.DM([[5, 7], [10, 3]]), ca.DM([[3, 4], [3, 2]]))])
    self.check("copysign", ca.Function("f", [x, y], [ca.copysign(x, y)]),
               [(ca.DM([[5, -3], [-7, 2]]), ca.DM([[-1, 1], [1, -1]]))])
    self.check("fmin", ca.Function("f", [x, y], [ca.fmin(x, y)]),
               [(ca.DM([[1, 5], [3, 2]]), ca.DM([[4, 2], [1, 6]]))])
    self.check("fmax", ca.Function("f", [x, y], [ca.fmax(x, y)]),
               [(ca.DM([[1, 5], [3, 2]]), ca.DM([[4, 2], [1, 6]]))])

  def test_linear_algebra(self):
    self.message("ONNX bilinear form and rank-1 update")
    A = ca.MX.sym("A", 3, 3)
    x = ca.MX.sym("x", 3, 1)
    y = ca.MX.sym("y", 3, 1)
    self.check("bilin", ca.Function("f", [A, x, y], [ca.bilin(A, x, y)]),
               [(ca.DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), ca.DM([1, 2, 3]), ca.DM([1, 0, 1]))])

    A = ca.MX.sym("A", 2, 2)
    alpha = ca.MX.sym("alpha")
    x = ca.MX.sym("x", 2, 1)
    y = ca.MX.sym("y", 2, 1)
    self.check("rank1", ca.Function("f", [A, alpha, x, y], [ca.rank1(A, alpha, x, y)]),
               [(ca.DM([[1, 2], [3, 4]]), ca.DM(2.0), ca.DM([1, 0]), ca.DM([0, 1]))])

  def test_atan2(self):
    self.message("ONNX atan2 across all quadrants and edge cases")
    y = ca.MX.sym("y")
    x = ca.MX.sym("x")
    self.check("atan2", ca.Function("f", [y, x], [ca.atan2(y, x)]),
               [(ca.DM(1.0), ca.DM(1.0)), (ca.DM(1.0), ca.DM(-1.0)),
                (ca.DM(-1.0), ca.DM(-1.0)), (ca.DM(-1.0), ca.DM(1.0)),
                (ca.DM(1.0), ca.DM(0.0)), (ca.DM(-1.0), ca.DM(0.0)),
                (ca.DM(0.0), ca.DM(1.0)), (ca.DM(0.0), ca.DM(-1.0)), (ca.DM(0.0), ca.DM(0.0))])

  def test_if_else(self):
    self.message("ONNX if_else_zero conditional")
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    self.check("if_else_zero", ca.Function("f", [x, y], [ca.if_else(x > 0, y, 0)]),
               [(ca.DM(1.0), ca.DM(5.0)), (ca.DM(-1.0), ca.DM(5.0)), (ca.DM(0.0), ca.DM(3.0))])

  def test_map(self):
    self.message("ONNX Map exported as Scan")
    x = ca.MX.sym("x")
    step = ca.Function("step", [x], [ca.sin(x) + 1])
    c = ca.MX.sym("c", 1, 3)
    self.check("map_scalar", ca.Function("f", [c], [step.map(3)(c)]), [ca.DM([[0.1, 0.5, 1.0]])])

    x = ca.MX.sym("x", 2)
    y = ca.MX.sym("y", 2)
    base = ca.Function("base", [x, y], [x + y, x * y])
    X = ca.MX.sym("X", 2, 3)
    Y = ca.MX.sym("Y", 2, 3)
    self.check("map_multi", ca.Function("f", [X, Y], base.map(3)(X, Y)),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM([[1, 1, 1], [2, 2, 2]]))])

  def test_map_multicolumn(self):
    self.message("ONNX Map over a multi-column base (3-D Scan envelope)")
    # base consumes a 2x2 block per iteration, produces a 2x1 column
    A = ca.MX.sym("A", 2, 2)
    f = ca.Function("f", [A], [ca.sum2(A)])
    X = ca.MX.sym("X", 2, 6)
    self.check("map_block_in", ca.Function("g", [X], [f.map(3)(X)]),
               [ca.DM([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]])])

    # base with multi-column input and output
    A = ca.MX.sym("A", 2, 2)
    f = ca.Function("f", [A], [A * 2])
    X = ca.MX.sym("X", 2, 4)
    self.check("map_block_inout", ca.Function("g", [X], [f.map(2)(X)]),
               [ca.DM([[1, 2, 3, 4], [5, 6, 7, 8]])])

  def test_map_reduce(self):
    self.message("ONNX reduce-map: captured (reduce_in) and accumulated (reduce_out) via Scan")
    # reduce_in: p is broadcast over the iterations
    x = ca.MX.sym("x", 2)
    p = ca.MX.sym("p", 1)
    base = ca.Function("base", [x, p], [x + p])
    X = ca.MX.sym("X", 2, 3)
    P = ca.MX.sym("P", 1)
    self.check("reduce_in", ca.Function("f", [X, P], [base.map(3, [False, True], [False])(X, P)]),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM(10))])

    # reduce_out: the per-iteration outputs are summed
    base = ca.Function("base", [x], [ca.sum1(x)])
    self.check("reduce_out", ca.Function("f", [X], [base.map(3, [False], [True])(X)]),
               [ca.DM([[1, 2, 3], [4, 5, 6]])])

    # both, with multiple outputs
    base = ca.Function("base", [x, p], [x + p, ca.sum1(x) * p])
    self.check("reduce_both", ca.Function("f", [X, P], list(base.map(3, [False, True], [False, True])(X, P))),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM(10))])

  def test_extra_math(self):
    self.message("ONNX det, logsumexp, log1p, expm1, hypot")
    A = ca.MX.sym("A", 2, 2)
    self.check("det", ca.Function("f", [A], [ca.det(A)]), [ca.DM([[1, 2], [3, 5]])])

    v = ca.MX.sym("v", 3, 1)
    self.check("logsumexp", ca.Function("f", [v], [ca.logsumexp(v)]), [ca.DM([0.1, 0.5, 1.0])])

    x = ca.MX.sym("x")
    self.check("log1p", ca.Function("f", [x], [ca.log1p(x)]), [ca.DM(0.5), ca.DM(2.0)])
    self.check("expm1", ca.Function("f", [x], [ca.expm1(x)]), [ca.DM(0.5), ca.DM(-1.0)])

    y = ca.MX.sym("y")
    self.check("hypot", ca.Function("f", [x, y], [ca.hypot(x, y)]), [(ca.DM(3.0), ca.DM(4.0))])

  def test_scatter(self):
    self.message("ONNX setnonzeros as ScatterElements")
    x = ca.MX.sym("x", 5)
    z = ca.MX.sym("z", 2)
    y = x + 0
    y[[1, 3]] = z
    self.check("scatter_indices", ca.Function("f", [x, z], [y]),
               [(ca.DM([1, 2, 3, 4, 5]), ca.DM([10, 20]))])

    y = x + 0
    y[1:3] = z
    self.check("scatter_slice", ca.Function("f", [x, z], [y]),
               [(ca.DM([1, 2, 3, 4, 5]), ca.DM([10, 20]))])

  def test_if(self):
    self.message("ONNX if_else (Switch) as If")
    x = ca.MX.sym("x", 2)
    ft = ca.Function("ft", [x], [x * 2])
    ff = ca.Function("ff", [x], [x + 1])
    sw = ca.Function.if_else("sw", ft, ff)
    c = ca.MX.sym("c")
    xx = ca.MX.sym("xx", 2)
    self.check("if_else", ca.Function("f", [c, xx], [sw(c, xx)]),
               [(ca.DM(1), ca.DM([3, 4])), (ca.DM(0), ca.DM([3, 4]))])

  def test_reshape(self):
    self.message("ONNX reshape (CasADi column-major vs ONNX row-major)")
    X = ca.MX.sym("X", 2, 3)
    self.check("reshape_3x2", ca.Function("g", [X], [ca.reshape(X, 3, 2)]), [ca.DM([[1, 2, 3], [4, 5, 6]])])
    self.check("reshape_6x1", ca.Function("g", [X], [ca.vec(X)]), [ca.DM([[1, 2, 3], [4, 5, 6]])])
    self.check("reshape_of_expr", ca.Function("g", [X], [ca.reshape(ca.sin(X), 3, 2)]),
               [ca.DM([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])])

  def test_index_2d(self):
    self.message("ONNX 2-D slicing / flat indexing (column-major)")
    M = ca.MX.sym("M", 3, 4)
    Mv = ca.DM([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    self.check("slice2d_rc", ca.Function("g", [M], [M[1:3, 1:4]]), [Mv])
    self.check("slice2d_cols", ca.Function("g", [M], [M[:, 1:3]]), [Mv])
    self.check("slice2d_rows", ca.Function("g", [M], [M[0:2, :]]), [Mv])
    self.check("index2d_flat", ca.Function("g", [M], [M[[0, 5, 10, 2]]]), [Mv])

  def test_scatter_2d(self):
    self.message("ONNX setnonzeros into a 2-D matrix")
    M = ca.MX.sym("M", 3, 4)
    Mv = ca.DM([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    z = ca.MX.sym("z", 2, 1)
    w = M + 0
    w[1:3, 1] = z
    self.check("scatter2d_col", ca.Function("g", [M, z], [w]), [(Mv, ca.DM([70, 80]))])
    zb = ca.MX.sym("zb", 2, 2)
    w = M + 0
    w[0:2, 1:3] = zb
    self.check("scatter2d_block", ca.Function("g", [M, zb], [w]), [(Mv, ca.DM([[70, 80], [90, 100]]))])

  def test_einstein(self):
    self.message("ONNX einstein contraction as Einsum (round-trip)")
    A = ca.MX.sym("A", 2, 3)
    B = ca.MX.sym("B", 3, 2)
    self.check("einstein_matmul",
               ca.Function("g", [A, B],
                           [ca.einstein(ca.vec(A), ca.vec(B), [2, 3], [3, 2], [2, 2],
                                        [-1, -3], [-3, -2], [-1, -2])]),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM([[1, 0], [0, 1], [1, 1]]))])

    A = ca.MX.sym("A", 2, 2)
    B = ca.MX.sym("B", 2, 2)
    self.check("einstein_inner",
               ca.Function("g", [A, B],
                           [ca.einstein(ca.vec(A), ca.vec(B), [2, 2], [2, 2], [],
                                        [-1, -2], [-1, -2], [])]),
               [(ca.DM([[1, 2], [3, 4]]), ca.DM([[5, 6], [7, 8]]))])

    u = ca.MX.sym("u", 3)
    v = ca.MX.sym("v", 2)
    self.check("einstein_outer",
               ca.Function("g", [u, v],
                           [ca.einstein(ca.vec(u), ca.vec(v), [3], [2], [3, 2],
                                        [-1], [-2], [-1, -2])]),
               [(ca.DM([1, 2, 3]), ca.DM([10, 20]))])

    A = ca.MX.sym("A", 2, 3)
    x = ca.MX.sym("x", 3)
    self.check("einstein_matvec",
               ca.Function("g", [A, x],
                           [ca.einstein(ca.vec(A), ca.vec(x), [2, 3], [3], [2],
                                        [-1, -2], [-2], [-1])]),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM([1, 1, 1]))])

  def test_unsupported_rejected(self):
    self.message("ONNX rejects multi-case Switch and mapaccum on export")
    x = ca.MX.sym("x")

    # Multi-case Switch (ONNX If is strictly 2-way)
    f0 = ca.Function("f0", [x], [x])
    f1 = ca.Function("f1", [x], [x * 2])
    fd = ca.Function("fd", [x], [x + 1])
    sw = ca.Function.conditional("sw", [f0, f1], fd)
    idx = ca.MX.sym("idx")
    xx = ca.MX.sym("xx")
    g = ca.Function("g", [idx, xx], [sw(idx, xx)])
    with self.assertRaises(RuntimeError):
      ca.GraphBuilder(g).export_onnx("unsupported.onnx")

    # mapaccum lowers to a multi-segment-I/O function (not yet supported as a nested function)
    u = ca.MX.sym("u")
    step = ca.Function("step", [x, u], [x + u])
    acc = step.mapaccum("acc", 3)
    x0 = ca.MX.sym("x0")
    U = ca.MX.sym("U", 1, 3)
    g = ca.Function("g", [x0, U], [acc(x0, U)])
    with self.assertRaises(RuntimeError):
      ca.GraphBuilder(g).export_onnx("unsupported.onnx")



if __name__ == '__main__':
  unittest.main()
