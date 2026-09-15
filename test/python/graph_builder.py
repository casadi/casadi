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
from typing import Union
from helpers import casadiTestCase, memory_heavy, requires_onnxbackend

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
  "norm_1", "norm_2", "norm_fro",   # ReduceL1/ReduceL2
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
@requires_onnxbackend("ort")
@unittest.skipUnless(have_onnx, "needs the onnx python package")
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
    inferred = ca.GraphBuilder(path).create("inferred", ["x", "p"], ["y"])
    self.assertEqual(inferred.is_diff_in(), [True, False])
    f = inferred
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

  def test_infer_diff_sibling_signatures(self):
    self.message("Forward and adjoint signatures infer input masks; explicit flags take precedence")
    path = self.isdiff_fwd_model(3)
    model = onnx.load(path)
    sibling = os.path.join(os.path.dirname(path), "fwd_" + os.path.basename(path))
    onnx.save(model, sibling)
    self._tmp.append(sibling)
    del model.graph.node[2:]
    del model.graph.input[2:]
    del model.graph.output[1:]
    onnx.save(model, path)
    f = ca.GraphBuilder(path).create("f")
    self.assertEqual(f.is_diff_in(), [True, False])
    self.assertEqual(ca.Function.deserialize(f.serialize()).is_diff_in(), [True, False])
    explicit = ca.GraphBuilder(path).create("explicit", {"is_diff_in": [True, True]})
    self.assertEqual(explicit.is_diff_in(), [True, True])
    with self.assertRaisesRegex(RuntimeError, "incomplete fwd signature"):
      explicit.forward(1)
    vi = lambda name: helper.make_tensor_value_info(name, TensorProto.FLOAT, [3, 1])
    graph = helper.make_graph([helper.make_node("Identity", ["adj_y"], ["adj_p"])],
      "conflict", [vi("adj_y")], [vi("adj_p")])
    adj = os.path.join(os.path.dirname(path), "adj_" + os.path.basename(path))
    onnx.save(helper.make_model(graph), adj)
    self._tmp.append(adj)
    with self.assertRaisesRegex(RuntimeError, "signatures disagree"):
      ca.GraphBuilder(path).create("conflict")
    os.remove(sibling)
    os.remove(adj)
    self.assertEqual(ca.GraphBuilder(path).create("primal").is_diff_in(), [True, True])

  def test_export_graph_model_filename(self):
    self.message("Graph call nodes identify their source ONNX files after serialization")
    import json
    path = self.affine_model(numpy.float32, (3,), 2, 1)
    sibling = self.sibling_model(path, "adj")
    f = ca.GraphBuilder(path).create("renamed")
    self.assertEqual(f.info()["model_path"], path)
    x = ca.MX.sym("x", 3)
    wrapper = ca.Function("wrapper", [x], [f(x), ca.jacobian(f(x), x)])
    for function in [wrapper, ca.Function.deserialize(wrapper.serialize())]:
      graph = json.loads(function.export_graph())
      calls = [node for part in [graph]+graph.get("functions", [])
               for node in part["nodes"] if "model_path" in node]
      self.assertEqual({os.path.normpath(node["model_path"]) for node in calls},
                       {os.path.normpath(path), os.path.normpath(sibling)})
      for node in calls:
        self.assertIn(os.path.basename(node["model_path"]), node["label"])

  def test_nondifferentiable_sibling_hessian(self):
    self.message("A runtime parameter needs no adjoint output or forward seed, including Hessians")
    def value(name, columns: Union[int, str]=1):
      return helper.make_tensor_value_info(name, TensorProto.FLOAT, [1, columns])
    with tempfile.TemporaryDirectory() as folder:
      def save(name, nodes, inputs, outputs, constants=()):
        graph = helper.make_graph(nodes, name, inputs, outputs, list(constants))
        model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
        model.ir_version = 8
        onnx.checker.check_model(model)
        onnx.save(model, os.path.join(folder, name+".onnx"))
      two = numpy_helper.from_array(numpy.array([[2]], dtype=numpy.float32), "two")
      save("f", [helper.make_node("Mul", ["x", "x"], ["xx"]),
        helper.make_node("Mul", ["p", "xx"], ["y"])], [value("x"), value("p")], [value("y")])
      save("adj_f", [helper.make_node("Mul", ["two", "p"], ["tp"]),
        helper.make_node("Mul", ["tp", "x"], ["tpx"]),
        helper.make_node("Mul", ["tpx", "adj_y"], ["adj_x"])],
        [value("x"), value("p"), value("adj_y", "nadj")], [value("adj_x", "nadj")], [two])
      flat = numpy_helper.from_array(numpy.array([1, -1], dtype=numpy.int64), "flat")
      save("fwd_adj_f", [
        helper.make_node("Transpose", ["fwd_x"], ["vx"], perm=[1, 0]),
        helper.make_node("MatMul", ["vx", "adj_y"], ["outer"]),
        helper.make_node("Reshape", ["outer", "flat"], ["packed"]),
        helper.make_node("Mul", ["x", "fwd_adj_y"], ["xdw"]),
        helper.make_node("Add", ["packed", "xdw"], ["sum"]),
        helper.make_node("Mul", ["two", "p"], ["tp"]),
        helper.make_node("Mul", ["tp", "sum"], ["fwd_adj_x"])],
        [value("x"), value("p"), value("adj_y", "nadj"), value("fwd_x", "nfwd"),
         value("fwd_adj_y", "nsens")], [value("fwd_adj_x", "nsens")], [two, flat])
      for symbolic in [False, True]:
        opts = {"symbolic": True, "is_diff_in": [True, False]} if symbolic else {}
        f = ca.GraphBuilder(os.path.join(folder, "f.onnx")).create("f", opts)
        self.assertEqual(f.is_diff_in(), [True, False])
        x, p = ca.MX.sym("x"), ca.MX.sym("p")
        y = f(x, p)
        h = ca.Function("h", [x, p], [ca.hessian(y, x)[0], ca.jacobian(y, p)])
        self.assertEqual(h.sparsity_out(1).nnz(), 0)
        for weather in [1, 3]:
          self.checkarray(2*weather, h(2, weather)[0])
          rev = f.reverse(2)
          self.assertEqual(rev.is_diff_in(), [True, False, True, True])
          self.assertEqual(rev.is_diff_out(), [True, False])
          adj = rev(2, weather, f(2, weather), ca.DM([[1, 2]]))
          self.checkarray(ca.DM([[4*weather, 8*weather]]), adj[0])
          self.checkarray(ca.DM.zeros(1, 2), adj[1])
        restored = ca.Function.deserialize(h.serialize())
        self.checkarray(6, restored(2, 3)[0])

  def sibling_model(self, path, kind, n=3):
    sources = {"fwd": self.fwd_model, "adj": self.rev_model, "jac": self.jac_model}
    model = onnx.load(sources[kind](n))
    graph = model.graph
    inputs = [t for t in graph.input if t.name.startswith(kind + "_")]
    outputs = [t for t in graph.output if t.name.startswith(kind + "_")]
    initializers = [t for t in graph.initializer if t.name != "two1"]
    graph.CopyFrom(helper.make_graph(list(graph.node)[1:], kind, inputs, outputs, initializers))
    sibling = os.path.join(os.path.dirname(path), kind + "_" + os.path.basename(path))
    onnx.checker.check_model(model)
    onnx.save(model, sibling)
    self._tmp.append(sibling)
    return sibling

  def test_sibling_derivatives(self):
    self.message("primal-only models discover dynamically seeded sibling derivatives")
    x = ca.DM([1, 2, 3])
    for kind in ["fwd", "adj", "jac"]:
      path = self.affine_model(numpy.float32, (3,), 2, 1)
      self.sibling_model(path, kind)
      f = ca.GraphBuilder(path).create("renamed")
      self.assertEqual(f.name_in(), ["x"])
      self.assertEqual(f.name_out(), ["y"])
      self.checkarray(2*x+1, f(x), "primal")
      # Serialize before constructing derivatives: sibling references must survive.
      f = ca.Function.deserialize(f.serialize())
      if kind == "jac":
        df = f.jacobian()
        self.assertEqual(df.sparsity_in("out_y"), ca.Sparsity(3, 1))
        self.checkarray(2*ca.DM.eye(3), df(x, 2*x+1), "sibling jacobian")
      else:
        for count in [1, 2, 8]:
          seeds = ca.reshape(ca.DM(list(range(3*count))), 3, count)
          df = f.forward(count) if kind == "fwd" else f.reverse(count)
          self.assertEqual(df.name(), kind + str(count) + "_renamed")
          self.assertEqual(df.sparsity_in("out_y"), ca.Sparsity(3, 1))
          df = ca.Function.deserialize(df.serialize())
          self.assertEqual(df.sparsity_in("out_y"), ca.Sparsity(3, 1))
          self.checkarray(2*seeds, df(x, 2*x+1, seeds), "sibling sensitivities")
      X = ca.MX.sym("X", 3)
      J = ca.Function("J", [X], [ca.jacobian(f(X), X)])
      self.checkarray(2*ca.DM.eye(3), J(x), "AD using sibling")
      self.assertNotIn("renamed", [g.name() for g in J.find_functions()])

  def test_sibling_uses_output(self):
    self.message("explicit nominal outputs retain their sparsity and reach the derivative")
    vi = lambda name, shape: helper.make_tensor_value_info(name, TensorProto.FLOAT, shape)
    for kind, seed, result, dim in [("fwd", "fwd_x", "fwd_y", "nfwd"),
                                   ("adj", "adj_y", "adj_x", "nadj")]:
      path = self.save_model(helper.make_graph(
          [helper.make_node("Exp", ["x"], ["y"])], "exp",
          [vi("x", [3])], [vi("y", [3])]))
      axes = numpy_helper.from_array(numpy.array([1], numpy.int64), "axes")
      sibling = self.save_model(helper.make_graph(
          [helper.make_node("Unsqueeze", ["out_y", "axes"], ["column"]),
           helper.make_node("Mul", ["column", seed], [result])], kind,
          [vi("x", [3]), vi("out_y", [3]), vi(seed, [3, dim])],
          [vi(result, [3, dim])], [axes]))
      destination = os.path.join(os.path.dirname(path), kind+"_"+os.path.basename(path))
      os.rename(sibling, destination)
      self._tmp.append(destination)
      f = ca.GraphBuilder(path).create("f")
      df = f.forward(2) if kind == "fwd" else f.reverse(2)
      self.assertEqual(df.sparsity_in("out_y"), ca.Sparsity.dense(3, 1))
      df = ca.Function.deserialize(df.serialize())
      seeds = ca.DM([[1, 2], [3, 4], [5, 6]])
      y = ca.DM([2, 3, 4])
      self.checkarray(ca.repmat(y, 1, 2)*seeds, df(ca.DM.zeros(3), y, seeds),
                      "supplied output is used")

  def test_sibling_nonlinear(self):
    self.message("sibling derivatives receive primal inputs for nonlinear sensitivities")
    vi = lambda name, shape: helper.make_tensor_value_info(name, TensorProto.FLOAT, shape)
    path = self.save_model(helper.make_graph(
        [helper.make_node("Mul", ["x", "x"], ["y"])], "square",
        [vi("x", [3])], [vi("y", [3])]))
    axes = numpy_helper.from_array(numpy.array([1], numpy.int64), "axes")
    two = numpy_helper.from_array(numpy.array(2, numpy.float32), "two")
    sibling = self.save_model(helper.make_graph(
        [helper.make_node("Unsqueeze", ["x", "axes"], ["column"]),
         helper.make_node("Mul", ["column", "two"], ["scale"]),
         helper.make_node("Mul", ["scale", "fwd_x"], ["fwd_y"])], "forward",
        [vi("x", [3]), vi("fwd_x", [3, "nfwd"])], [vi("fwd_y", [3, "nfwd"])],
        [axes, two]))
    destination = os.path.join(os.path.dirname(path), "fwd_" + os.path.basename(path))
    os.rename(sibling, destination)
    self._tmp.append(destination)
    f = ca.GraphBuilder(path).create("f")
    x = ca.DM([1, -2, 3])
    seeds = ca.DM([[1, 2], [3, 4], [5, 6]])
    self.checkarray(ca.repmat(2*x, 1, 2)*seeds, f.forward(2)(x, x*x, seeds), "nonlinear forward")
    X = ca.MX.sym("X", 3)
    J = ca.Function("J", [X], [ca.jacobian(f(X), X)])
    self.checkarray(ca.diag(2*x), J(x), "nonlinear AD")

  def test_sibling_lazy_loading(self):
    self.message("explicit masks keep sibling loading lazy; constructed derivatives embed their model")
    path = self.affine_model(numpy.float32, (3,), 2, 1)
    sibling = self.sibling_model(path, "fwd")
    with open(sibling, "wb") as fh:
      fh.write(b"invalid ONNX")
    with self.assertRaisesRegex(RuntimeError, "Failed to parse ONNX"):
      ca.GraphBuilder(path).create("infer")
    f = ca.GraphBuilder(path).create("f", {"is_diff_in": [True]})
    self.checkarray(ca.DM([3, 5, 7]), f(ca.DM([1, 2, 3])), "no derivative parsing")
    with self.assertRaises(RuntimeError):
      f.forward(2)
    self.sibling_model(path, "fwd")
    df = f.forward(2)
    serialized = df.serialize()
    os.remove(sibling)
    os.remove(path)
    df = ca.Function.deserialize(serialized)
    seeds = ca.DM([[1, 2], [3, 4], [5, 6]])
    self.checkarray(2*seeds, df(ca.DM.zeros(3), ca.DM.zeros(3), seeds), "embedded derivative")

  def test_sibling_lookup(self):
    self.message("sibling lookup uses the source path and ignores count-specific filenames")
    path = self.affine_model(numpy.float32, (3,), 2, 1)
    sibling = self.sibling_model(path, "fwd")
    numbered = os.path.join(os.path.dirname(path), "fwd8_" + os.path.basename(path))
    os.rename(sibling, numbered)
    self._tmp.append(numbered)
    f = ca.GraphBuilder(path).create("f")
    with self.assertRaises(RuntimeError):
      f.forward(8)
    os.rename(numbered, sibling)
    cwd = os.getcwd()
    try:
      os.chdir(os.path.dirname(path))
      b = ca.GraphBuilder(os.path.basename(path))
      os.chdir(cwd)
      f = b.create("different_name")
      df = f.forward(2)
      self.checkarray(2*ca.DM.ones(3, 2), df(ca.DM.zeros(3), ca.DM.zeros(3), ca.DM.ones(3, 2)),
                      "relative source after changing directory")
    finally:
      os.chdir(cwd)

  def test_sibling_symlink(self):
    self.message("sibling lookup follows the supplied filename when the primal is a symlink")
    if os.name == "nt":
      self.skipTest("creating symlinks may require administrator privileges")
    path = self.affine_model(numpy.float32, (3,), 2, 1)
    alias = path + ".alias.onnx"
    os.symlink(path, alias)
    self._tmp.append(alias)
    self.sibling_model(alias, "fwd")
    f = ca.GraphBuilder(alias).create("f")
    seeds = ca.DM.ones(3, 2)
    self.checkarray(2*seeds, f.forward(2)(ca.DM.zeros(3), ca.DM.zeros(3), seeds), "alias sibling")

  def test_sibling_validation(self):
    self.message("embedded derivatives take precedence; invalid sibling signatures fail clearly")
    path = self.fwd_model(3)
    sibling = self.sibling_model(path, "fwd")
    with open(sibling, "wb") as fh:
      fh.write(b"invalid ONNX")
    f = ca.GraphBuilder(path).create("f", ["x"], ["y"])
    self.checkarray(2*ca.DM.ones(3, 2), f.forward(2)(ca.DM.zeros(3), ca.DM.zeros(3), ca.DM.ones(3, 2)),
                    "embedded derivative precedence")
    for failure in ["signature", "shape"]:
      path = self.affine_model(numpy.float32, (3,), 2, 1)
      sibling = self.sibling_model(path, "fwd", n=4 if failure == "shape" else 3)
      if failure == "signature":
        model = onnx.load(sibling)
        model.graph.output[0].name = "wrong"
        model.graph.node[0].output[0] = "wrong"
        onnx.save(model, sibling)
      f = ca.GraphBuilder(path).create("f")
      with self.assertRaisesRegex(RuntimeError, "incomplete fwd signature|incompatible shape"):
        f.forward(2)

  def hessian_models(self, embedded=False):
    # f(x) = x0**3 + x1**3 + x0*x1, with a nonconstant, nondiagonal Hessian.
    vi = lambda name, shape: helper.make_tensor_value_info(name, TensorProto.DOUBLE, shape)
    node = helper.make_node
    constants = [numpy_helper.from_array(numpy.array([[0, 1], [1, 0]], numpy.float64), "A"),
                 numpy_helper.from_array(numpy.array(0.5), "half"),
                 numpy_helper.from_array(numpy.array(3.0), "three"),
                 numpy_helper.from_array(numpy.array(6.0), "six"),
                 numpy_helper.from_array(numpy.array([0], numpy.int64), "axis0"),
                 numpy_helper.from_array(numpy.array([1], numpy.int64), "axis1"),
                 numpy_helper.from_array(numpy.array([2], numpy.int64), "axis2"),
                 numpy_helper.from_array(numpy.array([2, -1], numpy.int64), "packed")]
    common = [node("Mul", ["x", "x"], ["square"]), node("MatMul", ["A", "x"], ["Ax"])]
    primal = common + [node("Mul", ["square", "x"], ["cube"]),
                       node("Mul", ["x", "Ax"], ["cross2"]),
                       node("Mul", ["cross2", "half"], ["cross"]),
                       node("Add", ["cube", "cross"], ["terms"]),
                       node("ReduceSum", ["terms", "axis0"], ["y"], keepdims=1)]
    grad = [node("Mul", ["three", "square"], ["quadratic"]),
            node("Add", ["quadratic", "Ax"], ["gradient"])]
    adjoint = [node("Mul", ["gradient", "adj_y"], ["adj_x"])]
    path = self.save_model(helper.make_graph(
        primal + (grad + adjoint if embedded else []), "f",
        [vi("x", [2, 1])] + ([vi("adj_y", [1, "nadj"])] if embedded else []),
        [vi("y", [1, 1])] + ([vi("adj_x", [2, "nadj"])] if embedded else []), constants))
    adj_path = os.path.join(os.path.dirname(path), "adj_" + os.path.basename(path))
    if not embedded:
      tmp = self.save_model(helper.make_graph(common + grad + adjoint, "adj_f",
          [vi("x", [2, 1]), vi("adj_y", [1, "nadj"])], [vi("adj_x", [2, "nadj"])], constants))
      os.rename(tmp, adj_path)
      self._tmp.append(adj_path)
    forward = common + grad + [
        node("Mul", ["six", "x"], ["diagonal"]),
        node("Mul", ["diagonal", "fwd_x"], ["diag_v"]),
        node("MatMul", ["A", "fwd_x"], ["A_v"]),
        node("Add", ["diag_v", "A_v"], ["Hv"]),
        node("Unsqueeze", ["Hv", "axis2"], ["Hv3"]),
        node("Unsqueeze", ["adj_y", "axis1"], ["w3"]),
        node("Mul", ["Hv3", "w3"], ["weighted3"]),
        node("Reshape", ["weighted3", "packed"], ["weighted"]),
        node("Mul", ["gradient", "fwd_adj_y"], ["seed_term"]),
        node("Add", ["weighted", "seed_term"], ["fwd_adj_x"])]
    tmp = self.save_model(helper.make_graph(forward, "fwd_adj_f",
        [vi("x", [2, 1]), vi("adj_y", [1, "nadj"]), vi("fwd_x", [2, "nfwd"]),
         vi("fwd_adj_y", [1, "nsens"])], [vi("fwd_adj_x", [2, "nsens"])], constants))
    forward_path = os.path.join(os.path.dirname(adj_path), "fwd_" + os.path.basename(adj_path))
    os.rename(tmp, forward_path)
    self._tmp.append(forward_path)
    return path, forward_path

  @memory_heavy()
  def test_sibling_hessian(self):
    self.message("forward-over-adjoint siblings supply exact Hessians and derivatives of seeds")
    X = ca.MX.sym("x", 2)
    expr = X[0]**3 + X[1]**3 + X[0]*X[1]
    reference = ca.Function("reference", [X], [expr])
    href = ca.Function("href", [X], [ca.hessian(expr, X)[0]])
    for embedded in [False, True]:
      path, forward_path = self.hessian_models(embedded)
      f = ca.GraphBuilder(path).create("renamed", ["x"], ["y"])
      f = ca.Function.deserialize(f.serialize())
      H, grad = ca.hessian(f(X), X)
      h = ca.Function("h", [X], [H])
      for x in [ca.DM([1.2, -0.7]), ca.DM([-0.3, 2.1])]:
        self.checkarray(href(x), h(x), "Hessian through forward-over-adjoint", digits=12)
        for nadj in [1, 2]:
          adj = f.reverse(nadj)
          adj_ref = reference.reverse(nadj)
          w = ca.DM([[1.5, -0.4]])[:, :nadj]
          adj_value = adj(x, f(x), w)
          for nfwd in [1, 3]:
            v = ca.reshape(ca.DM(list(range(1, 2*nfwd+1))), 2, nfwd)
            dw = ca.DM(numpy.arange(nadj*nfwd).reshape(1, -1) + 0.25)
            inputs = [x, f(x), w, adj_value, v, ca.DM.zeros(1, nfwd), dw]
            self.checkarray(adj_ref.forward(nfwd)(*inputs), adj.forward(nfwd)(*inputs),
                            "compounded sensitivities including adjoint seed", digits=12)
      # A constructed Hessian stays evaluable after its model files disappear.
      serialized = h.serialize()
      for source in list(self._tmp):
        if os.path.exists(source):
          os.remove(source)
      self.checkarray(href(ca.DM([1.2, -0.7])), ca.Function.deserialize(serialized)(ca.DM([1.2, -0.7])),
                      "serialized Hessian", digits=12)

  @memory_heavy()
  def test_sibling_reverse_over_reverse(self):
    self.message("adj_adj siblings keep independent seed dimensions at each reverse level")
    path, forward_path = self.hessian_models()
    os.remove(forward_path)  # Hessians must use reverse-over-reverse exclusively.
    vi = lambda name, shape: helper.make_tensor_value_info(name, TensorProto.DOUBLE, shape)
    node = helper.make_node
    constants = [numpy_helper.from_array(numpy.array([[0, 1], [1, 0]], numpy.float64), "A"),
                 numpy_helper.from_array(numpy.array(3.0), "three"),
                 numpy_helper.from_array(numpy.array(6.0), "six"),
                 numpy_helper.from_array(numpy.array([0], numpy.int64), "axis0"),
                 numpy_helper.from_array(numpy.array([2], numpy.int64), "axis2")]
    # Unpack outer reverse seeds as [input, outer direction, inner direction].
    # Infer the outer count from the packed seed width and inner count through Reshape.
    constants.append(numpy_helper.from_array(numpy.array([2, -1], numpy.int64), "head"))
    constants.append(numpy_helper.from_array(numpy.array([1], numpy.int64), "axis1"))
    nodes = [node("Shape", ["adj_y"], ["wshape"]),
             node("Gather", ["wshape", "axis1"], ["inner"]),
             node("Concat", ["head", "inner"], ["shape3"], axis=0),
             node("Reshape", ["adj2_adj_x", "shape3"], ["u3"]),
             node("Unsqueeze", ["adj_y", "axis2"], ["wcol"]),
             node("MatMul", ["u3", "wcol"], ["weighted3"]),
             node("Squeeze", ["weighted3", "axis2"], ["weighted"]),
             node("Mul", ["six", "x"], ["diagonal"]),
             node("Mul", ["diagonal", "weighted"], ["diag_term"]),
             node("MatMul", ["A", "weighted"], ["offdiag_term"]),
             node("Add", ["diag_term", "offdiag_term"], ["adj2_x"]),
             node("Mul", ["x", "x"], ["square"]),
             node("Mul", ["three", "square"], ["quadratic"]),
             node("MatMul", ["A", "x"], ["Ax"]),
             node("Add", ["quadratic", "Ax"], ["gradient"]),
             node("Mul", ["gradient", "adj2_adj_x"], ["products"]),
             node("ReduceSum", ["products", "axis0"], ["adj2_adj_y"], keepdims=1)]
    tmp = self.save_model(helper.make_graph(nodes, "adj_adj_f",
        [vi("x", [2, 1]), vi("adj_y", [1, "nadj"]),
         vi("adj2_adj_x", [2, "packed_reverse"])],
        [vi("adj2_x", [2, "nadj2"]), vi("adj2_adj_y", [1, "packed_reverse"])], constants))
    sibling = os.path.join(os.path.dirname(path), "adj_adj_" + os.path.basename(path))
    os.rename(tmp, sibling)
    self._tmp.append(sibling)
    X = ca.MX.sym("x", 2)
    expr = X[0]**3 + X[1]**3 + X[0]*X[1]
    reference = ca.Function("reference", [X], [expr])
    f = ca.Function.deserialize(ca.GraphBuilder(path).create("f").serialize())
    for ninner in [1, 2, 3]:
      for nouter in [1, 2, 4]:
        a, ar = f.reverse(ninner), reference.reverse(ninner)
        x = ca.DM([1.2, -0.7])
        w = ca.DM(numpy.arange(ninner).reshape(1, -1) + 0.3)
        u = ca.DM(numpy.arange(2*ninner*nouter).reshape(2, -1) - 0.7)
        args = [x, f(x), w, a(x, f(x), w), u]
        actual = ca.Function.deserialize(a.reverse(nouter).serialize())(*args)
        for expected, result in zip(ar.reverse(nouter)(*args), actual):
          self.checkarray(expected, result, "reverse-over-reverse", digits=11)
    h = ca.Function("h", [X], [ca.hessian(f(X), X)[0]])
    href = ca.Function("href", [X], [ca.hessian(expr, X)[0]])
    self.checkarray(href(x), h(x), "Hessian using only reverse siblings", digits=12)

  @memory_heavy()
  def test_sibling_deep_nesting(self):
    self.message("three successive forward or reverse derivatives use recursive sibling lookup")
    for mode in ["adj", "fwd"]:
      X = ca.MX.sym("x")
      reference = ca.Function("reference", [X], [X**4], ["x"], ["y"])
      with tempfile.TemporaryDirectory() as directory:
        filename = "f.onnx"
        refs = [reference]
        for depth in range(4):
          target = os.path.join(directory, filename)
          ca.GraphBuilder(refs[-1]).export_onnx(target)
          # Omit structural-zero derivative placeholders from the ONNX interface.
          model = onnx.load(target)
          for io, count, nnz in [(model.graph.input, refs[-1].n_in(), refs[-1].nnz_in),
                                 (model.graph.output, refs[-1].n_out(), refs[-1].nnz_out)]:
            for i in reversed(range(count)):
              if nnz(i) == 0:
                del io[i]
          onnx.save(model, target)
          if depth < 3:
            refs.append(refs[-1].reverse(1) if mode == "adj" else refs[-1].forward(1))
          filename = mode + "_" + filename
        actual = ca.GraphBuilder(os.path.join(directory, "f.onnx")).create("f")
        for depth, ref in enumerate(refs):
          if depth:
            actual = actual.reverse(1) if mode == "adj" else actual.forward(1)
          actual = ca.Function.deserialize(actual.serialize())
          args = [ca.DM.ones(ref.sparsity_in(i))*((i+1)*0.3) for i in range(ref.n_in())]
          for expected, result in zip(ref.call(args), actual.call(args)):
            self.checkarray(expected, result, "nested " + mode, digits=12)

  def test_sibling_lazy_compound_lookup(self):
    self.message("compounded derivative files are looked up on demand without a cached file list")
    path, forward_path = self.hessian_models()
    delayed = forward_path + ".delayed"
    os.rename(forward_path, delayed)
    self._tmp.append(delayed)
    f = ca.GraphBuilder(path).create("f")
    serialized = f.serialize()
    os.rename(delayed, forward_path)
    f = ca.Function.deserialize(serialized)
    X = ca.MX.sym("x", 2)
    H = ca.Function("H", [X], [ca.hessian(f(X), X)[0]])
    self.checkarray(ca.DM([[6, 1], [1, 12]]), H(ca.DM([1, 2])), "late Hessian model")

  @memory_heavy()
  def test_sibling_configuration(self):
    self.message("nested siblings inherit specialized dimensions, baked inputs and seed names")
    vi = lambda name, shape: helper.make_tensor_value_info(name, TensorProto.DOUBLE, shape)
    node = helper.make_node
    # y = p*x**2 gives a nonzero Hessian after specializing batch and baking p.
    path = self.save_model(helper.make_graph(
        [node("Mul", ["x", "x"], ["square"]), node("Mul", ["square", "p"], ["y"])],
        "primal", [vi("x", ["batch", 1]), vi("p", [])], [vi("y", ["batch", 1])]))
    constants = [numpy_helper.from_array(numpy.array(2.0), "two"),
                 numpy_helper.from_array(numpy.array([1], numpy.int64), "axis1"),
                 numpy_helper.from_array(numpy.array([2], numpy.int64), "axis2"),
                 numpy_helper.from_array(numpy.array([0, -1], numpy.int64), "packed")]
    gradient = [node("Mul", ["two", "p"], ["scale"]),
                node("Mul", ["scale", "x"], ["gradient"])]
    for kind, seed, result, dim in [("fwd", "fwd_x", "fwd_y", "directions"),
                                    ("adj", "adj_y", "adj_x", "weights")]:
      sibling = self.save_model(helper.make_graph(
          gradient + [node("Mul", ["gradient", seed], [result])], kind,
          [vi("x", ["batch", 1]), vi("p", []), vi(seed, ["batch", dim])],
          [vi(result, ["batch", dim])], constants))
      destination = os.path.join(os.path.dirname(path), kind + "_" + os.path.basename(path))
      os.rename(sibling, destination)
      self._tmp.append(destination)
    mixed = gradient + [
        node("Mul", ["scale", "fwd_x"], ["Hv"]),
        node("Unsqueeze", ["Hv", "axis2"], ["Hv3"]),
        node("Unsqueeze", ["adj_y", "axis1"], ["w3"]),
        node("Mul", ["Hv3", "w3"], ["weighted3"]),
        node("Reshape", ["weighted3", "packed"], ["weighted"]),
        node("Mul", ["gradient", "fwd_adj_y"], ["seed_term"]),
        node("Add", ["weighted", "seed_term"], ["fwd_adj_x"])]
    sibling = self.save_model(helper.make_graph(mixed, "fwd_adj",
        [vi("x", ["batch", 1]), vi("p", []), vi("adj_y", ["batch", "weights"]),
         vi("fwd_x", ["batch", "directions"]), vi("fwd_adj_y", ["batch", "nsens"])],
        [vi("fwd_adj_x", ["batch", "nsens"])], constants))
    destination = os.path.join(os.path.dirname(path), "fwd_adj_" + os.path.basename(path))
    os.rename(sibling, destination)
    self._tmp.append(destination)
    for batch in [2, 3]:
      b = ca.GraphBuilder(path)
      b.bind_dim("batch", batch)
      b.set("p", 5)
      f = b.create("f", {"fwd_dim": "directions", "adj_dim": "weights"})
      del b
      f = ca.Function.deserialize(f.serialize())
      X = ca.MX.sym("x", batch)
      reference = ca.Function("reference", [X], [5*X**2])
      h = ca.Function("h", [X], [ca.hessian(ca.sum1(f(X)), X)[0]])
      for x in [ca.DM(range(1, batch+1)), -ca.DM(range(2, batch+2))]:
        self.checkarray(reference(x), f(x), "specialized primal", digits=12)
        self.checkarray(10*ca.DM.eye(batch), h(x), "specialized Hessian", digits=12)
        for nfwd in [1, 3]:
          v = ca.reshape(ca.DM(range(1, batch*nfwd+1)), batch, nfwd)
          self.checkarray(reference.forward(nfwd)(x, reference(x), v),
                          f.forward(nfwd)(x, f(x), v), "configured forward", digits=12)
          for nadj in [1, 2]:
            w = ca.reshape(ca.DM(range(1, batch*nadj+1)), batch, nadj)
            adj = f.reverse(nadj)
            adj_ref = reference.reverse(nadj)
            value = adj(x, f(x), w)
            self.checkarray(adj_ref(x, reference(x), w), value, "configured reverse", digits=12)
            dw = ca.reshape(ca.DM(range(1, batch*nadj*nfwd+1)), batch, nadj*nfwd)
            inputs = [x, f(x), w, value, v, ca.DM.zeros(batch, nfwd), dw]
            self.checkarray(adj_ref.forward(nfwd)(*inputs), adj.forward(nfwd)(*inputs),
                            "configured forward-over-reverse including seed derivatives", digits=12)

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

  def test_finite_difference_masks(self):
    self.message("FD preserves runtime parameters and nondifferentiable outputs")
    vi = lambda name: helper.make_tensor_value_info(name, TensorProto.DOUBLE, [1])
    path = self.save_model(helper.make_graph([
        helper.make_node("Mul", ["x", "x"], ["xx"]),
        helper.make_node("Mul", ["p", "xx"], ["y"]),
        helper.make_node("Add", ["x", "p"], ["z"])], "masked_fd",
        [vi("x"), vi("p")], [vi("y"), vi("z")]))
    f = ca.GraphBuilder(path).create("f", {"enable_fd": True, "fd_options": {"h": 1e-3},
        "is_diff_in": [True, False], "is_diff_out": [True, False]})
    self.assertEqual(f.forward(1).class_name(), "CentralDiff")
    x = ca.MX.sym("x")
    p = ca.MX.sym("p")
    y, z = f(x, p)
    D = ca.Function("D", [x, p], [y, z,
        ca.jacobian(ca.vertcat(y, z), ca.vertcat(x, p)), ca.hessian(y, x)[0]])
    for restored in [D, ca.Function.deserialize(D.serialize())]:
      for value in [1., 3.]:
        result = restored(2, value)
        self.checkarray(result[0], 4*value)
        self.checkarray(result[1], 2+value)
        self.checkarray(result[2], ca.DM([[4*value, 0], [0, 0]]), digits=4)
        self.checkarray(result[3], 2*value, digits=3)

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
    self.assertEqual(list(b.dimension("x")), [2, 3])
    self.assertEqual(b.dtype("x"), "FLOAT")
    self.assertEqual(list(b.dimension_param("x")), ["", ""])
    self.assertEqual(list(b.dynamic_params()), [])

  def test_builder_dynamic_dim(self):
    self.message("GraphBuilder dynamic-dimension binding")
    path = self.dyn_model()
    b = ca.GraphBuilder(path)
    self.assertEqual(list(b.dynamic_params()), ["batch"])
    self.assertEqual(list(b.dimension("x")), [-1, 3])
    self.assertEqual(list(b.dimension_param("x")), ["batch", ""])
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
    self.assertIn('#include <onnxruntime_c_api.h>', code)

  def test_codegen_compile(self):
    # Full roundtrip: generate C, compile+link against ONNX Runtime, load and evaluate.
    # Gated on ONNXRUNTIME_DIR (path to an ONNX Runtime C/C++ distribution).
    ort_dir = os.environ.get("ONNXRUNTIME_DIR")
    if not ort_dir:
      self.skipTest("set ONNXRUNTIME_DIR to a C/C++ distribution to run the codegen compile test")
    import subprocess
    self.message("generated C compiles against ONNX Runtime and evaluates correctly")
    path = self.affine_model(numpy.float64, (3,), 2.0, 1.0)
    sibling = self.sibling_model(path, "fwd")
    f = ca.GraphBuilder(path).create("f")
    df = f.forward(2)
    os.remove(sibling)
    x = ca.DM([0.5, 1.0, -2.0])
    seeds = ca.DM([[1, 2], [3, 4], [5, 6]])
    hessian_path, _ = self.hessian_models()
    hf = ca.GraphBuilder(hessian_path).create("hf")
    hx = ca.MX.sym("x", 2)
    h = ca.Function("h", [hx], [ca.hessian(hf(hx), hx)[0]])
    cases = [(f, [x], 2*x+1), (df, [x, 2*x+1, seeds], 2*seeds),
             (h, [ca.DM([1, 2])], ca.DM([[6, 1], [1, 12]]))]
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
      try:
        os.chdir(d)
        for fn, inputs, expected in cases:
          fn.generate(fn.name() + ".c")
          inc = os.path.join(os.path.dirname(ca.__file__), "include", "casadi", "interfaces", "ort")
          r = subprocess.run(["gcc", "-shared", "-fPIC", fn.name() + ".c",
                              "-I" + inc, "-I" + os.path.join(ort_dir, "include"),
                              "-L" + os.path.join(ort_dir, "lib"), "-lonnxruntime",
                              "-o", fn.name() + ".so"], capture_output=True, text=True)
          self.assertEqual(r.returncode, 0, r.stderr)
          g = ca.external(fn.name(), os.path.join(d, fn.name() + ".so"))
          self.checkarray(expected, g(*inputs), "codegen eval", digits=10)
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

  def test_builder_print(self):
    self.message("GraphBuilder printing")
    b = ca.GraphBuilder(self.two_pathway_model())
    self.assertTrue("GraphBuilder" in str(b) and "input(s)" in str(b))

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

  @unittest.skipUnless(have_onnx, "needs the onnx python package")
  def test_constant_numeric_attributes(self):
    self.message("ONNX scalar and vector Constant attributes become doubles")
    cases = [
      ("value_int", -7, TensorProto.INT64, [], [-7]),
      ("value_int", 2**40, TensorProto.INT64, [], [2**40]),
      ("value_ints", [-3, 0, 7], TensorProto.INT64, [3], [-3, 0, 7]),
      ("value_float", 1.25, TensorProto.FLOAT, [], [1.25]),
      ("value_floats", [-2.5, 0., 1.25], TensorProto.FLOAT, [3], [-2.5, 0., 1.25]),
    ]
    with tempfile.TemporaryDirectory() as folder:
      path = os.path.join(folder, "constant.onnx")
      for attr, value, dtype, shape, expected in cases:
        with self.subTest(attribute=attr, value=value):
          node = helper.make_node("Constant", [], ["y"], **{attr: value})
          graph = helper.make_graph([node], "constant", [],
            [helper.make_tensor_value_info("y", dtype, shape)])
          model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
          onnx.checker.check_model(model)
          onnx.save(model, path)
          f = ca.GraphBuilder(path).create("f", {"symbolic": True})
          self.checkarray(ca.DM(expected), f()["y"], attr)

  @unittest.skipUnless(have_onnx, "needs the onnx python package")
  def test_constant_attribute_reshape_ad(self):
    self.message("Integer Constant attributes drive symbolic reshape and AD")
    nodes = [
      helper.make_node("Constant", [], ["shape"], value_ints=[2, 3]),
      helper.make_node("Reshape", ["x", "shape"], ["r"]),
      helper.make_node("Constant", [], ["scale"], value_int=3),
      helper.make_node("Cast", ["scale"], ["scale_double"], to=TensorProto.DOUBLE),
      helper.make_node("Mul", ["r", "scale_double"], ["y"]),
    ]
    graph = helper.make_graph(nodes, "reshape", [
      helper.make_tensor_value_info("x", TensorProto.DOUBLE, [6])], [
      helper.make_tensor_value_info("y", TensorProto.DOUBLE, [2, 3])])
    model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
    onnx.checker.check_model(model)
    with tempfile.TemporaryDirectory() as folder:
      path = os.path.join(folder, "reshape.onnx")
      onnx.save(model, path)
      f = ca.GraphBuilder(path).create("f", {"symbolic": True})
      x = ca.MX.sym("x", f.sparsity_in(0))
      ref = ca.Function("ref", [x], [3*ca.reshape(x, 3, 2).T])
      self.checkfunction(ref, f, inputs=[ca.DM([1, 2, 3, 4, 5, 6])])

  @unittest.skipUnless(have_onnx, "needs the onnx python package")
  def test_constant_slice_bounds(self):
    self.message("ONNX constant Slice bounds use regular getnonzeros")
    cases = [
      (1, 2**63-1, 1, [1, 2, 3, 4, 5]),
      (-2**63, 4, 1, [0, 1, 2, 3]),
      (1, 2**63-1, 2, [1, 3, 5]),
      (2**63-1, -2**63, -1, [5, 4, 3, 2, 1, 0]),
      (-2, -2**63, -2, [4, 2, 0]),
      (1, 2**63-1, 2**63-1, [1]),
      (5, -2**63, -2**63, [5]),
      (4, 1, 1, []),
      (-2**63, -2**63, -1, []),
    ]
    with tempfile.TemporaryDirectory() as folder:
      path = os.path.join(folder, "slice.onnx")
      for start, end, step, expected_indices in cases:
        with self.subTest(start=start, end=end, step=step):
          nodes = [helper.make_node("Constant", [], [name], value_ints=[value])
            for name, value in [("start", start), ("end", end), ("axis", 0), ("step", step)]]
          nodes.append(helper.make_node("Slice", ["x", "start", "end", "axis", "step"], ["y"]))
          graph = helper.make_graph(nodes, "slice", [
            helper.make_tensor_value_info("x", TensorProto.DOUBLE, [6, 1])], [
            helper.make_tensor_value_info("y", TensorProto.DOUBLE, [len(expected_indices), 1])])
          model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
          onnx.checker.check_model(model)
          onnx.save(model, path)
          f = ca.GraphBuilder(path).create("f", {"symbolic": True})
          self.checkarray(ca.DM(expected_indices), f(ca.DM(range(6))))
          ops = [f.instruction_id(k) for k in range(f.n_instructions())]
          self.assertNotIn(ca.OP_GETNONZEROS_PARAM, ops)
          if expected_indices:
            self.assertIn(ca.OP_GETNONZEROS, ops)

  @unittest.skipUnless(have_onnx, "needs the onnx python package")
  def test_external_primal_symbolic_ad(self):
    self.message("External ONNX nominal network preserves shapes and uses CasADi AD")
    w = numpy.array([[1., 2.], [-1., .5], [.25, -.75]])
    b = numpy.array([.1, -.2, .3])
    nodes = [helper.make_node("Constant", [], [name], value_ints=value)
      for name, value in [("start", [0]), ("end", [2**63-1]), ("axis", [0]),
        ("shape", [1, 2]), ("flat", [-1]), ("unsqueeze_axis", [1])]]
    nodes += [
      helper.make_node("Constant", [], ["index"], value_int=0),
      helper.make_node("Slice", ["x", "start", "end", "axis"], ["sliced"]),
      helper.make_node("Gather", ["sliced", "index"], ["selected"], axis=1),
      helper.make_node("Reshape", ["selected", "shape"], ["row"]),
      helper.make_node("Gemm", ["row", "w", "b"], ["affine"]),
      helper.make_node("Tanh", ["affine"], ["activation"]),
      helper.make_node("Reshape", ["activation", "flat"], ["vector"]),
      helper.make_node("Unsqueeze", ["vector", "unsqueeze_axis"], ["y"]),
    ]
    graph = helper.make_graph(nodes, "network", [
      helper.make_tensor_value_info("x", TensorProto.DOUBLE, [2, 1])], [
      helper.make_tensor_value_info("y", TensorProto.DOUBLE, [3, 1])], [
      numpy_helper.from_array(w.T.copy(), "w"), numpy_helper.from_array(b, "b")])
    model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
    onnx.checker.check_model(model)
    with tempfile.TemporaryDirectory() as folder:
      path = os.path.join(folder, "f.onnx")
      onnx.save(model, path)
      f = ca.GraphBuilder(path).create("f", {"symbolic": True})
      self.assertEqual(f.size_in(0), (2, 1))
      self.assertEqual(f.size_out(0), (3, 1))
      x = ca.MX.sym("x", 2)
      ref = ca.Function("ref", [x], [ca.tanh(ca.mtimes(w, x) + b)])
      self.checkfunction(ref, f, inputs=[ca.DM([.2, -.1])])
      h = ca.Function("h", [x], [ca.hessian(ca.sumsqr(f(x)), x)[0]])
      href = ca.Function("href", [x], [ca.hessian(ca.sumsqr(ref(x)), x)[0]])
      self.checkarray(href([.2, -.1]), h([.2, -.1]))
      ops = [f.instruction_id(k) for k in range(f.n_instructions())]
      self.assertNotIn(ca.OP_GETNONZEROS_PARAM, ops)

  @unittest.skipUnless(have_onnx, "needs the onnx python package")
  def test_external_gather_constant_indices(self):
    self.message("External ONNX Gather uses fixed indices along the declared axis")
    values = numpy.arange(6.).reshape(2, 3)
    with tempfile.TemporaryDirectory() as folder:
      path = os.path.join(folder, "gather.onnx")
      for axis in [0, 1]:
        for indices in [1, [1, 0]]:
          expected = numpy.take(values, indices, axis=axis)
          if isinstance(indices, int):
            constant = helper.make_node("Constant", [], ["indices"], value_int=indices)
          else:
            constant = helper.make_node("Constant", [], ["indices"], value_ints=indices)
          nodes = [constant, helper.make_node("Gather", ["x", "indices"], ["y"], axis=axis)]
          graph = helper.make_graph(nodes, "gather", [
            helper.make_tensor_value_info("x", TensorProto.DOUBLE, [2, 3])], [
            helper.make_tensor_value_info("y", TensorProto.DOUBLE, list(expected.shape))])
          model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
          onnx.checker.check_model(model)
          onnx.save(model, path)
          f = ca.GraphBuilder(path).create("f", {"symbolic": True})
          self.assertEqual(f.size_in(0), (2, 3))
          self.checkarray(ca.DM(expected), f(values))
          ops = [f.instruction_id(k) for k in range(f.n_instructions())]
          self.assertIn(ca.OP_GETNONZEROS, ops)
          self.assertNotIn(ca.OP_GETNONZEROS_PARAM, ops)

  @unittest.skipUnless(have_onnx, "needs the onnx python package")
  def test_direct_integer_structural_inputs(self):
    self.message("Direct integer constants feed Slice, Reshape and Gather without rounding")
    with tempfile.TemporaryDirectory() as folder:
      path = os.path.join(folder, "integers.onnx")
      for dtype in [numpy.int32, numpy.int64]:
        for encoding in ["attribute", "tensor", "raw", "initializer", "raw_initializer"]:
          with self.subTest(dtype=dtype, encoding=encoding):
            nodes, initializers = [], []
            def constant(name, values, scalar=False):
              shape = [] if scalar else [len(values)]
              if encoding == "attribute":
                nodes.append(helper.make_node("Constant", [], [name], **{
                  "value_int" if scalar else "value_ints": values[0] if scalar else values}))
                return
              if "raw" in encoding:
                tensor = numpy_helper.from_array(numpy.array(values, dtype=dtype).reshape(shape), name)
              else:
                elem_type = TensorProto.INT32 if dtype == numpy.int32 else TensorProto.INT64
                tensor = helper.make_tensor(name, elem_type, shape, values)
              if "initializer" in encoding:
                initializers.append(tensor)
              else:
                nodes.append(helper.make_node("Constant", [], [name], value=tensor))
            constant("start", [-6])
            constant("end", [numpy.iinfo(dtype).max])
            constant("axis", [0])
            constant("step", [2])
            constant("shape", [1, 3])
            constant("index", [-2], scalar=True)
            nodes += [
              helper.make_node("Slice", ["x", "start", "end", "axis", "step"], ["s"]),
              helper.make_node("Reshape", ["s", "shape"], ["r"]),
              helper.make_node("Gather", ["r", "index"], ["y"], axis=1),
            ]
            graph = helper.make_graph(nodes, "integers", [
              helper.make_tensor_value_info("x", TensorProto.DOUBLE, [6, 1])], [
              helper.make_tensor_value_info("y", TensorProto.DOUBLE, [1])], initializers)
            model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
            onnx.checker.check_model(model)
            onnx.save(model, path)
            f = ca.GraphBuilder(path).create("f", {"symbolic": True})
            self.checkarray(ca.DM(2), f(ca.DM(range(6))))
            ops = [f.instruction_id(k) for k in range(f.n_instructions())]
            self.assertIn(ca.OP_GETNONZEROS, ops)
            self.assertNotIn(ca.OP_GETNONZEROS_PARAM, ops)
            if dtype != numpy.int64:
              continue
            # An invalid large index/shape must be reported exactly, not rounded through double.
            for op in ["Gather", "Reshape"]:
              nodes, initializers = [], []
              exact = 2**53 + 1
              constant("integer", [exact], scalar=op == "Gather")
              nodes.append(helper.make_node(op, ["x", "integer"], ["y"]))
              graph = helper.make_graph(nodes, "exact", [
                helper.make_tensor_value_info("x", TensorProto.DOUBLE, [1])], [
                helper.make_tensor_value_info("y", TensorProto.DOUBLE, [1])], initializers)
              model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
              onnx.checker.check_model(model)
              onnx.save(model, path)
              with self.assertRaisesRegex(RuntimeError, str(exact)):
                ca.GraphBuilder(path).create("f", {"symbolic": True})

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
    # Transpose-representation: every ONNX tensor is the transpose of its CasADi value, so feed
    # the transposed inputs and un-transpose ORT's outputs before comparing.
    fn = "test_%s_ort.onnx" % name
    dtype = numpy.float32 if float32 else numpy.float64
    try:
      ca.GraphBuilder(f).export_onnx(fn, {"casadi_real": "float" if float32 else "double"})
      session = ort.InferenceSession(fn)
      input_names = [i.name for i in session.get_inputs()]
      for vals in inputs:
        if not isinstance(vals, tuple):
          vals = (vals,)
        feed = {input_names[i]: numpy.array(ca.DM(v)).T.astype(dtype) for i, v in enumerate(vals)}
        ref = f(*vals)
        if not isinstance(ref, (list, tuple)):
          ref = [ref]
        for a, b in zip(ref, session.run(None, feed)):
          self.checkarray(a, ca.DM(numpy.array(b).T), name, digits=digits)
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

  def test_sparse_constant(self):
    # A sparse DM constant must survive CasADi -> ONNX -> CasADi as an ONNX sparse_value (COO),
    # preserving the sparsity PATTERN (not densifying it) and the values. The symbolic (import)
    # round-trip checks pattern preservation; ONNX Runtime also runs the model -- it densifies the
    # sparse Constant at its Constant node and yields the correct dense result, which we verify too.
    self.message("ONNX sparse constant (sparse_value / COO) round-trip")
    x = ca.MX.sym("x")
    cases = {
      "diag": ca.DM(ca.diag(ca.DM([2, 4, 5]))),                       # symmetric pattern
      "asym": ca.DM.triplet([0, 1, 2], [1, 0, 2], [7, 5, 3], 3, 3),   # asymmetric: stresses row/col order
      "rect": ca.DM.triplet([0, 2], [1, 0], [9, -3], 3, 4),           # non-square, off-diagonal
    }
    for name, A in cases.items():
      self.assertTrue(not A.is_dense(), "%s should be sparse" % name)
      f = ca.Function("f", [x], [A * x])
      fn = "test_sparse_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        # Stored faithfully as a sparse Constant (sparse_value). Transpose-representation: the
        # ONNX tensor is A^T, so dims are reversed and the COO indices are A^T's (col,row), still
        # ascending row-major and unique.
        if have_onnx:
          m = onnx.load(fn)
          sparse_consts = [a.sparse_tensor for n in m.graph.node if n.op_type == "Constant"
                           for a in n.attribute if a.name == "sparse_value"]
          self.assertEqual(len(sparse_consts), 1, "%s: expected one sparse Constant" % name)
          st = sparse_consts[0]
          self.assertEqual(list(st.dims), [A.size2(), A.size1()])
          idx = [tuple(p) for p in numpy_helper.to_array(st.indices).tolist()]
          self.assertEqual(len(idx), A.nnz())
          self.assertEqual(idx, sorted(idx), "%s: COO indices must be ascending (row,col)" % name)
          self.assertEqual(len(set(idx)), len(idx), "%s: COO indices must be unique" % name)
        g = ca.GraphBuilder(fn).create("imported_sparse_%s" % name, {"symbolic": True})
        # Pattern preserved (not densified) ...
        self.assertTrue(f.sparsity_out(0) == g.sparsity_out(0),
                        "%s: sparsity not preserved" % name)
        self.assertEqual(g.nnz_out(0), A.nnz())
        # ... and numerically identical
        for v in [2.0, -1.5]:
          self.checkarray(f(v), g(v), "sparse_const_%s" % name, digits=12)
        # ONNX Runtime executes the model (densifying the sparse Constant) -> correct dense result.
        # Feed the scalar input as rank-2 [[v]]: CasADi declares inputs 2-D, so a rank-0 scalar
        # would trip ORT's "Invalid rank for input ... Got: 0 Expected: 2" (unrelated to sparsity).
        if has_onnxruntime:
          session = ort.InferenceSession(fn)
          iname = session.get_inputs()[0].name
          for v in [2.0, -1.5]:
            out = session.run(None, {iname: numpy.array([[v]], dtype=numpy.float64)})[0]
            # ORT yields the transpose-representation; un-transpose before comparing
            self.checkarray(f(v), ca.DM(numpy.array(out).T), "sparse_const_ort_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_sparse_slice(self):
    # Slicing a SPARSE value (a single row / column) is an OP_GETNONZEROS whose output is itself
    # sparse. The exporter recognizes it as a structured row/col slice and emits a clean ONNX Slice
    # (not a flat Gather); the importer rebuilds A[r,:]/A[:,c] on the sparse value, so the SPARSITY
    # PATTERN round-trips -- no densify, no overlay needed.
    self.message("ONNX sparse row/column slice (getnonzeros -> Slice) round-trip")
    # intricate 4x4 pattern, distinct values so r stays an OP_CONST
    r = ca.DM.triplet([0, 1, 2, 2, 3, 1], [0, 0, 1, 3, 2, 3], [5, 7, 9, 11, 13, 15], 4, 4)
    x = ca.MX.sym("x")
    cases = {"row2": (r * x)[2, :], "row1": (r * x)[1, :],
             "col3": (r * x)[:, 3], "col0": (r * x)[:, 0]}
    for name, expr in cases.items():
      self.assertTrue(not expr.sparsity().is_dense(), "%s should be sparse" % name)
      f = ca.Function("f", [x], [expr])
      fn = "test_sparse_slice_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        if have_onnx:  # recognized as a structured slice, not a flat gather
          ops = [n.op_type for n in onnx.load(fn).graph.node]
          self.assertTrue("Slice" in ops, "%s: expected a Slice node" % name)
          # CLEAN: no flatten -- the sparse value flows straight through Mul -> Slice
          self.assertTrue("Gather" not in ops and "Reshape" not in ops,
                          "%s: should have no flatten (Gather/Reshape): %s" % (name, ops))
        g = ca.GraphBuilder(fn).create("imported_sparse_slice_%s" % name, {"symbolic": True})
        self.assertTrue(f.sparsity_out(0) == g.sparsity_out(0),
                        "%s: sparsity pattern not preserved" % name)
        self.assertEqual(g.nnz_out(0), expr.nnz())
        for v in [2.0, -1.5]:
          self.checkarray(f(v), g(v), "sparse_slice_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_dense_slice_recognizer(self):
    # On a DENSE input, an OP_GETNONZEROS that forms a regular grid A[r0:r1:rs, c0:c1:cs] -- single
    # row/col, range, block, or strided -- is recognized geometrically and emitted as a clean ONNX
    # Slice (the canonical subtensor op), not a flat Gather. Irregular index sets fall back to Gather.
    self.message("ONNX dense getnonzeros -> Slice recognizer (row/col/range/block/strided)")
    A = ca.MX.sym("A", 4, 5)
    av = ca.DM(numpy.arange(1, 21).reshape(5, 4).T)
    grid = {"row": A[2, :], "col": A[:, 3], "row_range": A[1:3, :], "col_range": A[:, 1:4],
            "block": A[1:3, 1:4], "strided_rows": A[0:4:2, :], "strided_cols": A[:, 0:5:2],
            "strided_block": A[1:3, 1:4:2]}
    for name, expr in grid.items():
      f = ca.Function("f", [A], [expr])
      fn = "test_slicerec_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        if have_onnx:
          ops = [n.op_type for n in onnx.load(fn).graph.node]
          self.assertTrue("Slice" in ops, "%s: expected a Slice node, got %s" % (name, ops))
          # CLEAN graph: a recognized slice is NOT flattened -- no Gather, no Reshape at all
          self.assertTrue("Gather" not in ops, "%s: should not flatten to Gather" % name)
          self.assertTrue("Reshape" not in ops, "%s: should have no flatten Reshape" % name)
        g = ca.GraphBuilder(fn).create("imported_slicerec_%s" % name, {"symbolic": True})
        self.checkarray(ca.DM(f(av)), ca.DM(g(av)), "slicerec_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)
    # an irregular index set is NOT a grid -> must fall back to the flat path (Reshape + Gather)
    firr = ca.Function("f", [A], [A.nz[[3, 1, 7]]])
    fn = "test_slicerec_irregular.onnx"
    try:
      ca.GraphBuilder(firr).export_onnx(fn)
      if have_onnx:
        ops = [n.op_type for n in onnx.load(fn).graph.node]
        self.assertTrue("Gather" in ops and "Slice" not in ops, "irregular should use Gather: %s" % ops)
        self.assertTrue("Reshape" in ops, "irregular fallback should flatten (Reshape): %s" % ops)
      g = ca.GraphBuilder(fn).create("imported_slicerec_irr", {"symbolic": True})
      self.checkarray(ca.DM(firr(av)), ca.DM(g(av)), "slicerec_irregular", digits=12)
    finally:
      if os.path.exists(fn):
        os.remove(fn)

  def test_getnonzeros_param(self):
    # y = x.nz[p] with p a RUNTIME index vector (OP_GETNONZEROS_PARAM, op 70) -- the parametric
    # twin of the constant getnonzeros. p[k] addresses x's column-major NONZERO array. Export
    # flattens x's pseudo-dense col-major to a 1 x numel row and GatherElements by the runtime
    # indices (mapping nz->dense col-major position via a constant find() table when x is sparse);
    # import reconstructs it with the MX get_nz overload. Output is a vector shaped like p.
    self.message("ONNX parametric getnonzeros x.nz[p] (OP_GETNONZEROS_PARAM) round-trip")
    x = ca.MX.sym("x", 5)
    p = ca.MX.sym("p", 3)
    A = ca.MX.sym("A", 3, 3)
    q = ca.MX.sym("q", 4)
    S = ca.MX.sym("S", ca.Sparsity.lower(3))
    r = ca.MX.sym("r", 2)
    cases = {
      "dense_vec": (ca.Function("f", [x, p], [x.nz[p]]),
                    (ca.DM(range(5)), ca.DM([0, 2, 4]))),
      "dense_mat": (ca.Function("f", [A, q], [A.nz[q]]),
                    (ca.DM(numpy.arange(9).reshape(3, 3)), ca.DM([0, 4, 8, 1]))),
      "sparse":    (ca.Function("f", [S, r], [S.nz[r]]),
                    (ca.DM(ca.Sparsity.lower(3), [1, 2, 3, 4, 5, 6]), ca.DM([0, 3]))),
    }
    for name, (f, vals) in cases.items():
      fn = "test_gnzp_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        if have_onnx:
          ops = [n.op_type for n in onnx.load(fn).graph.node]
          self.assertTrue("GatherElements" in ops,
                          "%s: expected a GatherElements node, got %s" % (name, ops))
        g = ca.GraphBuilder(fn).create("imported_gnzp_%s" % name, {"symbolic": True})
        self.assertIn(ca.OP_GETNONZEROS_PARAM,
          [g.instruction_id(k) for k in range(g.n_instructions())])
        self.assertTrue(f.sparsity_out(0) == g.sparsity_out(0),
                        "%s: sparsity pattern not preserved" % name)
        self.checkarray(ca.DM(f(*vals)), ca.DM(g(*vals)), "gnzp_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_setnonzeros_scatternd(self):
    # The write-dual of the Slice recognizer: setnonzeros (y.nz[idx]=z) emits a 2-D ONNX ScatterND
    # that writes z's values at the (row,col) COORDINATES of the target cells -- no flatten/reshape
    # of x. Each coordinate's column-major linear position recovers the original nz position, so the
    # import reconstructs it via set_nz on the 2-D value. Works for row/col/block/scalar/sparse-z.
    self.message("ONNX setnonzeros -> ScatterND (2-D coordinate write)")
    x = ca.MX.sym("x", 4, 5)
    xv = ca.DM(numpy.arange(20).reshape(5, 4).T)
    z_row, z_col = ca.MX.sym("zr", 1, 5), ca.MX.sym("zc", 4, 1)
    z_blk, z_sca = ca.MX.sym("zb", 2, 3), ca.MX.sym("zs")
    z_sp = ca.MX.sym("zsp", ca.Sparsity.diag(2))
    def mk(setter):
      Y = ca.MX(x + 0); setter(Y); return Y
    cases = {
      "row":    (mk(lambda Y: Y.__setitem__((1, slice(None)), z_row)), z_row, ca.DM([[1, 2, 3, 4, 5]])),
      "col":    (mk(lambda Y: Y.__setitem__((slice(None), 2), z_col)), z_col, ca.DM([[6], [7], [8], [9]])),
      "block":  (mk(lambda Y: Y.__setitem__((slice(1, 3), slice(1, 4)), z_blk)), z_blk, ca.DM([[1, 2, 3], [4, 5, 6]])),
      "scalar": (mk(lambda Y: Y.__setitem__((2, 3), z_sca)), z_sca, ca.DM(99)),
      "sparse_z": (mk(lambda Y: Y.__setitem__((slice(1, 3), slice(1, 3)), z_sp)), z_sp,
                   ca.DM(ca.Sparsity.diag(2), [7, 9])),
    }
    for name, (Y, zsym, zval) in cases.items():
      f = ca.Function("f", [x, zsym], [Y])
      fn = "test_scnd_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        if have_onnx:
          m = onnx.load(fn)
          ops = [n.op_type for n in m.graph.node]
          self.assertTrue("ScatterND" in ops, "%s: expected ScatterND, got %s" % (name, ops))
          # CLEAN graph: the scatter TARGET (x) is NOT flattened -- it flows straight into ScatterND.
          producer = {o: n.op_type for n in m.graph.node for o in n.output}
          sc = [n for n in m.graph.node if n.op_type == "ScatterND"][0]
          self.assertTrue(producer.get(sc.input[0]) != "Reshape",
                          "%s: ScatterND target x was flattened (Reshape)" % name)
          # dense-z writes use ONLY ScatterND (no flat ScatterElements). The sparse_z case has an
          # incidental ScatterElements from the internal getnonzeros the sparse-block assign induces.
          if name != "sparse_z":
            self.assertTrue("ScatterElements" not in ops, "%s: should not use flat ScatterElements" % name)
        g = ca.GraphBuilder(fn).create("imported_scnd_%s" % name, {"symbolic": True})
        self.checkarray(ca.DM(f(xv, zval)), ca.DM(g(xv, zval)), "scnd_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_sparse_fancy_index(self):
    # getnonzeros cases the row/col recognizer does NOT catch (a sparse row RANGE, a sparse BLOCK)
    # fall through to the general fallback, which must NOT crash: it scatters the gathered nonzeros
    # into a dense zero frame. The ONNX intermediate is dense, but the sparsity overlay restores
    # the output pattern on import (see test_sparsity_overlay); numerics round-trip exactly.
    self.message("ONNX sparse fancy indexing fallback (range/block) numeric round-trip")
    r = ca.DM.triplet([0, 1, 2, 2, 3, 1], [0, 0, 1, 3, 2, 3], [5, 7, 9, 11, 13, 15], 4, 4)
    x = ca.MX.sym("x")
    for name, expr in {"row_range": (r * x)[1:3, :], "block": (r * x)[1:3, 1:3],
                       "col_range": (r * x)[:, 2:4]}.items():
      self.assertTrue(not expr.sparsity().is_dense(), "%s should be sparse" % name)
      f = ca.Function("f", [x], [expr])
      fn = "test_sparse_fancy_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)  # must not raise (sparse-output getnonzeros)
        g = ca.GraphBuilder(fn).create("imported_fancy_%s" % name, {"symbolic": True})
        for v in [2.0, -1.5]:
          self.checkarray(f(v), g(v), "sparse_fancy_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_project_densify(self):
    # OP_PROJECT / densify (a nonzero rearrangement between two sparsities) routes through the same
    # poorest-form remap helper as getnonzeros -- it must export (no "unsupported op 75") and
    # round-trip numerically. The ONNX intermediate is dense; the sparsity overlay restores a
    # non-dense output pattern on import (see test_sparsity_overlay).
    self.message("ONNX project / densify round-trip")
    A = ca.MX.sym("A", ca.Sparsity.diag(3))
    cases = {"densify": ca.densify(A),
             "project_dense": ca.project(A, ca.Sparsity.dense(3, 3)),
             "project_lower": ca.project(A, ca.Sparsity.lower(3))}
    vals = [ca.DM(ca.Sparsity.diag(3), [2, 4, 5]), ca.DM(ca.Sparsity.diag(3), [-1, 0.5, 3])]
    for name, expr in cases.items():
      f = ca.Function("f", [A], [expr])
      fn = "test_project_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)  # must not raise
        g = ca.GraphBuilder(fn).create("imported_project_%s" % name, {"symbolic": True})
        for v in vals:
          self.checkarray(f(v), g(v), "project_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_project_mask_downstream(self):
    # OP_PROJECT exports as a single Mul by a DENSE 0/1 mask, keeping the body uniformly
    # pseudo-dense (no sparse intermediate). A project result must still feed ANY downstream op
    # correctly -- slice, getnonzeros, setnonzeros, matmul, reproject -- with the I/O sparsity
    # signature restored by the overlay. Regression for project composed with downstream ops.
    self.message("ONNX project mask feeding downstream ops")
    B = ca.MX.sym("B", 3, 3)
    P = ca.project(B, ca.Sparsity.lower(3))           # sparse (lower-tri) intermediate
    z = ca.MX.sym("z", 2, 2)
    Yset = ca.MX(P); Yset[1:3, 0:2] = z               # project result feeds a setnonzeros
    cases = {
      "as_output":   ca.Function("f", [B], [P]),
      "row_slice":   ca.Function("f", [B], [P[2, :]]),
      "getnonzeros": ca.Function("f", [B], [P.nz[[1, 3, 5]]]),
      "matmul":      ca.Function("f", [B], [ca.mtimes(B, P)]),
      "setnz_block": ca.Function("f", [B, z], [Yset]),
      "reproject":   ca.Function("f", [B], [ca.project(P, ca.Sparsity.diag(3))]),
    }
    v = ca.DM(numpy.arange(1, 10).reshape(3, 3).T)
    zv = ca.DM([[7, 0], [0, 9]])
    for name, f in cases.items():
      args = [v, zv] if f.n_in() == 2 else [v]
      fn = "test_projdown_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        g = ca.GraphBuilder(fn).create("imported_projdown_%s" % name, {"symbolic": True})
        for o in range(f.n_out()):
          self.assertTrue(f.sparsity_out(o) == g.sparsity_out(o), "%s: out sparsity" % name)
        self.checkarray(ca.DM(f(*args)), ca.DM(g(*args)), "projdown_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_sparse_setnonzeros(self):
    # setnonzeros (y.nz[idx] = z) with a SPARSE update source z, and into a SPARSE target. The
    # write positions are z's nonzeros, so the exporter must feed z's NONZEROS as the updates
    # (matching the index count) -- else ScatterElements indices/updates mismatch on import. The
    # output densifies (pattern preservation is a follow-up) but numerics round-trip exactly.
    self.message("ONNX setnonzeros with sparse update source / sparse target")
    x = ca.MX.sym("x", 4, 4)
    zc = ca.MX.sym("zc", ca.Sparsity.diag(2))   # sparse RHS
    Yb = ca.MX(x + 0); Yb[1:3, 1:3] = zc
    self.message("setnz dense[1:3,1:3]=sparse2x2")
    self._roundtrip_numeric("setnz_block_sparse_z", ca.Function("f", [x, zc], [Yb]),
                            [(numpy.arange(16).reshape(4, 4), ca.DM(ca.Sparsity.diag(2), [7, 9]))])
    ssp = ca.MX.sym("ssp", ca.Sparsity.lower(3))
    zz = ca.MX.sym("zz")
    Ts = ca.MX(ssp); Ts[2, 0] = zz              # write an existing nonzero of a sparse target
    self._roundtrip_numeric("setnz_sparse_target", ca.Function("f", [ssp, zz], [Ts]),
                            [(ca.DM(ca.Sparsity.lower(3), [1, 2, 3, 4, 5, 6]), 9.0)])

  def test_addnonzeros(self):
    # ADDNONZEROS (z.nz[idx] += v, accumulating DUPLICATE idx) exports as ScatterND reduction="add"
    # and imports as the reverse-mode adjoint of a getnonzeros (it cannot be built as an MX directly).
    # An addnonzeros arises from reverse AD; build the cases with jtimes(getnonzeros, ..., tr=True).
    self.message("ONNX addnonzeros (scatter-add) via getnonzeros adjoint round-trip")
    x = ca.MX.sym("x", 5); s = ca.MX.sym("s", 4)
    A = ca.MX.sym("A", 3, 3); sm = ca.MX.sym("sm", 4)
    S = ca.MX.sym("S", ca.Sparsity.lower(3)); ss = ca.MX.sym("ss", 3)
    cases = {
      # repeated index 2 -> the two seeds accumulate; into a zero frame and onto a nonzero base
      "dupidx":    (ca.Function("f", [x, s], [ca.jtimes(x.nz[[0, 2, 4, 2]], x, s, True)]),
                    [(ca.DM(range(5)), ca.DM([1, 2, 3, 4]))]),
      "onto_base": (ca.Function("f", [x, s], [x + ca.jtimes(x.nz[[0, 2, 4, 2]], x, s, True)]),
                    [(ca.DM(range(5)), ca.DM([1, 2, 3, 4]))]),
      "matrix":    (ca.Function("f", [A, sm], [ca.jtimes(A.nz[[0, 4, 8, 4]], A, sm, True)]),
                    [(ca.DM(numpy.arange(9).reshape(3, 3)), ca.DM([1, 2, 3, 4]))]),
      "sparse":    (ca.Function("f", [S, ss], [ca.jtimes(S.nz[[0, 2, 5]], S, ss, True)]),
                    [(ca.DM(ca.Sparsity.lower(3), [1, 2, 3, 4, 5, 6]), ca.DM([7, 8, 9]))]),
    }
    for name, (f, inputs) in cases.items():
      fn = "test_addnz_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        g = ca.GraphBuilder(fn).create("imported_addnz_%s" % name, {"symbolic": True})
        self.assertTrue(f.sparsity_out(0) == g.sparsity_out(0), "%s: sparsity not preserved" % name)
        for vals in inputs:
          self.checkarray(ca.DM(f(*vals)), ca.DM(g(*vals)), "addnz_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_setnonzeros_param(self):
    # SETNONZEROS_PARAM (y.nz[p]=v) and ADDNONZEROS_PARAM (y.nz[p]+=v, accumulating DUPLICATE p) with
    # a RUNTIME index vector p. Export as ScatterElements on the flattened pseudo-dense target
    # (reduction "none"/"add"); import via the set_nz MX overload / the getnonzeros-adjoint (jtimes).
    # The add-param variant has no direct MX builder -> build it with jtimes of a getnonzeros_param.
    self.message("ONNX set/addnonzeros with a parametric (runtime) index vector")
    x = ca.MX.sym("x", 5); v = ca.MX.sym("v", 3); p = ca.MX.sym("p", 3); s = ca.MX.sym("s", 3)
    A = ca.MX.sym("A", 3, 3); vm = ca.MX.sym("vm", 2); pm = ca.MX.sym("pm", 2)
    S = ca.MX.sym("S", ca.Sparsity.lower(3)); vs = ca.MX.sym("vs", 2); ps = ca.MX.sym("ps", 2)
    yv = ca.MX(x); yv.nz[p] = v
    ym = ca.MX(A); ym.nz[pm] = vm
    yS = ca.MX(S); yS.nz[ps] = vs
    cases = {
      "set_vec":    (ca.Function("f", [x, v, p], [yv]),
                     [(ca.DM(range(5)), ca.DM([7, 8, 9]), ca.DM([0, 2, 4]))]),
      "set_mat":    (ca.Function("f", [A, vm, pm], [ym]),
                     [(ca.DM(numpy.arange(9).reshape(3, 3)), ca.DM([10, 20]), ca.DM([0, 8]))]),
      "set_sparse": (ca.Function("f", [S, vs, ps], [yS]),
                     [(ca.DM(ca.Sparsity.lower(3), [1, 2, 3, 4, 5, 6]), ca.DM([70, 80]), ca.DM([0, 5]))]),
      # add-param, runtime p with a repeated index -> the two seeds accumulate
      "add_dup":    (ca.Function("f", [x, p, s], [ca.jtimes(x.nz[p], x, s, True)]),
                     [(ca.DM(range(5)), ca.DM([0, 2, 2]), ca.DM([7, 8, 9]))]),
      "add_onbase": (ca.Function("f", [x, p, s], [x + ca.jtimes(x.nz[p], x, s, True)]),
                     [(ca.DM(range(5)), ca.DM([1, 3, 3]), ca.DM([7, 8, 9]))]),
    }
    for name, (f, inputs) in cases.items():
      fn = "test_snzp_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        g = ca.GraphBuilder(fn).create("imported_snzp_%s" % name, {"symbolic": True})
        self.assertTrue(f.sparsity_out(0) == g.sparsity_out(0), "%s: sparsity not preserved" % name)
        for vals in inputs:
          self.checkarray(ca.DM(f(*vals)), ca.DM(g(*vals)), "snzp_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_sparse_showcase(self):
    # End-to-end pseudo-dense showcase: one graph with several SPARSE inputs/constants and a mix of
    # ops (sparse+sparse union, sparse*dense mtimes, a sparse row slice, a pattern-changing project,
    # a runtime-index gather, a block-diagonal concat) producing six outputs with SIX DISTINCT
    # sparsity patterns. Asserts the contract on BOTH sides: the ONNX is dense value-flow with
    # sparsity only in SEEDS (sparse_value constants + a sparse_initializer, NO metadata overlay) and
    # clean op promotions (Gemm/Slice/GatherElements, no flatten-gather jungle); the imported MX
    # reproduces every output's exact pattern AND values (symbolically and through ONNX Runtime).
    self.message("ONNX sparse showcase: dense graph + seeds -> exact sparsity round-trip")
    A = ca.MX.sym("A", ca.Sparsity.lower(4))            # sparse symbolic input (10 nz)
    x = ca.MX.sym("x", 4); p = ca.MX.sym("p", 3)        # dense vector, runtime index vector
    D = ca.DM(ca.diag(ca.DM([10, 20, 30, 40])))         # sparse CONSTANT seed (4 nz)
    M = A + D                                           # sparse+sparse -> lower union (10 nz)
    expr = [M, ca.mtimes(M, x), M[2, :], ca.project(M, ca.Sparsity.diag(4)),
            (ca.mtimes(M, x)).nz[p], ca.diagcat(M, M[2, :])]
    names = ["M", "y", "row2", "diagM", "picked", "blk"]
    f = ca.Function("showcase", [A, x, p], expr, ["A", "x", "p"], names)
    # the test must be non-trivial: outputs span distinct, genuinely-sparse patterns
    want = {"M": (4, 4, 10), "y": (4, 1, 4), "row2": (1, 4, 3),
            "diagM": (4, 4, 4), "picked": (3, 1, 3), "blk": (5, 8, 13)}
    for i, nm in enumerate(names):
      s = f.sparsity_out(i)
      self.assertEqual((s.size1(), s.size2(), s.nnz()), want[nm], "%s shape/nnz" % nm)

    fn = "test_showcase.onnx"
    try:
      ca.GraphBuilder(f).export_onnx(fn)

      if have_onnx:
        mdl = onnx.load(fn)
        # opset 16 (ScatterND/-Elements reduction available)
        self.assertEqual({o.domain: o.version for o in mdl.opset_import}.get(""), 16)
        # sparsity lives ONLY in seeds: a sparse_initializer for the sparse input + sparse_value
        # constants; and emphatically NO metadata overlay (the deleted cruft).
        self.assertEqual(len(mdl.graph.sparse_initializer), 1, "one sparse input seed")
        n_sparse_const = sum(1 for n in mdl.graph.node if n.op_type == "Constant"
                             for a in n.attribute if a.name == "sparse_value")
        self.assertGreaterEqual(n_sparse_const, 1, "sparse_value constant seeds present")
        self.assertFalse(any("sparsity" in kv.key for kv in mdl.metadata_props),
                         "no casadi_sparsity_* metadata overlay")
        # clean op promotions, and value-flow is plain dense ONNX ops (no jungle / no custom ops)
        ops = set(n.op_type for n in mdl.graph.node)
        for need in ("Gemm", "Slice", "GatherElements"):  # mtimes / row-slice / runtime gather
          self.assertIn(need, ops, "expected clean %s" % need)
        allowed = {"Add", "Sub", "Mul", "Div", "Cast", "Constant", "Identity", "Gemm",
                   "Slice", "GatherElements", "ScatterElements", "Pad", "Sum", "Concat",
                   "Transpose", "Reshape", "Tile"}
        self.assertTrue(ops <= allowed, "unexpected (jungle?) ops: %s" % (ops - allowed))
        self.assertLessEqual(len(mdl.graph.node), 30, "graph stays terse (no jungle)")

      g = ca.GraphBuilder(fn).create("imported_showcase", {"symbolic": True})
      # MX side: identical I/O signature, every output's EXACT pattern preserved
      self.assertEqual((g.n_in(), g.n_out()), (f.n_in(), f.n_out()))
      for i in range(f.n_in()):
        self.assertTrue(f.sparsity_in(i) == g.sparsity_in(i), "input %d sparsity" % i)
      for i in range(f.n_out()):
        self.assertTrue(f.sparsity_out(i) == g.sparsity_out(i), "%s sparsity" % names[i])

      Av = ca.DM(ca.Sparsity.lower(4), list(range(1, 11)))
      xv = ca.DM([1, 2, 3, 4]); pv = ca.DM([0, 3, 1])
      ref, got = f(A=Av, x=xv, p=pv), g(A=Av, x=xv, p=pv)
      for nm in names:                                   # symbolic round-trip values
        self.checkarray(ref[nm], got[nm], "showcase_%s" % nm, digits=12)

      if has_onnxruntime:                                # and through ONNX Runtime (dense execution)
        sess = ort.InferenceSession(fn)
        feed = {"A": numpy.array(Av).T, "x": numpy.array(ca.DM(xv)).T.reshape(1, 4),
                "p": numpy.array(ca.DM(pv)).T.reshape(1, 3)}
        out = {o.name: v for o, v in zip(sess.get_outputs(), sess.run(None, feed))}
        for nm in names:
          self.checkarray(ref[nm], ca.DM(numpy.array(out[nm]).T), "showcase_ort_%s" % nm, digits=9)
    finally:
      if os.path.exists(fn):
        os.remove(fn)

  def test_passthrough_ops(self):
    # OP_ASSERTION (attachAssert) and OP_MONITOR are value-passthroughs (the assert condition /
    # monitor print are side effects); they export as Identity and must round-trip the data unchanged.
    self.message("ONNX assertion / monitor passthrough round-trip")
    x = ca.MX.sym("x", 3)
    cases = {
      "assertion": ca.Function("f", [x], [x.attachAssert(ca.dot(x, x) < 100, "m") * 2]),
      "monitor":   ca.Function("f", [x], [ca.sin(x.monitor("tag"))]),
    }
    for name, f in cases.items():
      fn = "test_pass_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        g = ca.GraphBuilder(fn).create("imported_pass_%s" % name, {"symbolic": True})
        self.assertTrue(f.sparsity_out(0) == g.sparsity_out(0), "%s: sparsity" % name)
        self.checkarray(ca.DM(f(ca.DM([1, 2, 3]))), ca.DM(g(ca.DM([1, 2, 3]))), name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_sparsity_overlay(self):
    # The ONNX graph is dense, but a non-dense CasADi I/O Sparsity round-trips via SEEDS, not
    # metadata: a sparse input is a sparse_initializer (sparse input symbol on import); a non-dense
    # result is restored in-graph by its producing op (mask Mul / zero Add), so CasADi recovers the
    # exact pattern natively even where an intermediate densifies. No metadata, no output overlay.
    self.message("ONNX sparse I/O signature round-trip (seed-based, no metadata)")
    r = ca.DM.triplet([0, 1, 2, 2, 3, 1], [0, 0, 1, 3, 2, 3], [5, 7, 9, 11, 13, 15], 4, 4)
    x = ca.MX.sym("x")
    xi = ca.MX.sym("xi", ca.Sparsity.diag(3))                  # sparse INPUT
    cases = {
      "sparse_input": (ca.Function("f", [xi], [xi * 2]), [ca.DM(ca.Sparsity.diag(3), [2, 4, 5])]),
      "fallback_row_range": (ca.Function("f", [x], [(r * x)[1:3, :]]), [3.0]),   # densified, seed-restored
      "project_lower": (ca.Function("f", [xi], [ca.project(xi, ca.Sparsity.lower(3))]),
                        [ca.DM(ca.Sparsity.diag(3), [2, 4, 5])]),
    }
    for name, (f, args) in cases.items():
      fn = "test_overlay_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        g = ca.GraphBuilder(fn).create("imported_overlay_%s" % name, {"symbolic": True})
        for i in range(f.n_in()):
          self.assertTrue(f.sparsity_in(i) == g.sparsity_in(i), "%s: input sparsity" % name)
        for o in range(f.n_out()):
          self.assertTrue(f.sparsity_out(o) == g.sparsity_out(o), "%s: output sparsity" % name)
        self.checkarray(ca.DM(f(*args)), ca.DM(g(*args)), name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_compound_sparse_chain(self):
    # End-to-end integration of the whole sparsity pipeline in one graph:
    #   lower-tri sparse INPUT -> square -> project onto diagonal -> setnonzeros of a row
    #   -> sum over rows -> sum over columns -> scalar.
    # Exercises sparse-input overlay, elementwise, project (Mul-mask), setnonzeros (ScatterND),
    # and reductions together. Plus the (r*x).nz[2,0,1] sparse irregular-permutation getnonzeros
    # (flat-Gather fallback on a sparse value).
    self.message("ONNX compound sparse chain (square/project/setnz/reductions)")
    x = ca.MX.sym("x", ca.Sparsity.lower(3))
    pd = ca.project(x ** 2, ca.Sparsity.diag(3))
    y = ca.MX(pd); y[1, :] = ca.DM([[10, 20, 30]])
    f = ca.Function("f", [x], [ca.sum2(ca.sum1(y))])
    # hand value: diag(x^2) = (2^2,5^2,7^2); row 1 -> [10,20,30]; sum = 4 + 49 + (10+20+30) = 113
    self._roundtrip_numeric("compound_chain", f,
                            [ca.DM(ca.Sparsity.lower(3), [2, 3, 4, 5, 6, 7])])
    self.assertEqual(float(f(ca.DM(ca.Sparsity.lower(3), [2, 3, 4, 5, 6, 7]))), 113.0)

    xs = ca.MX.sym("xs")
    r = ca.DM(ca.diag([2.0, 4.0, 5.0]))            # sparse constant (distinct values)
    fperm = ca.Function("f", [xs], [(r * xs).nz[[2, 0, 1]]])   # irregular nz permutation -> column
    self._roundtrip_numeric("nz_permutation", fperm, [3.0])

  def test_sparse_stress_edges(self):
    # Edge cases found by stress-testing the sparsity support. Notably an all-structural-zero
    # output: the MXFunction has NO instructions, so the exporter must emit an explicit zero
    # Constant for the output (else neither ORT nor the importer have a tensor to produce).
    self.message("ONNX sparsity stress edge cases")
    s = ca.MX.sym("s")
    S = ca.MX.sym("S", ca.Sparsity.lower(4))
    sv = ca.DM(ca.Sparsity.lower(4), range(1, 11))
    A = ca.MX.sym("A", 4, 5)
    av = ca.DM(numpy.arange(20).reshape(5, 4).T)
    cases = {
      "zero_output_0nnz": (ca.Function("f", [s], [ca.MX(ca.Sparsity(3, 3)) * s]), [2.0]),
      "1x1_empty_plus_scalar": (ca.Function("f", [s], [ca.MX(ca.Sparsity(1, 1)) + s]), [5.0]),
      "pattern_collapse_S_minus_S": (ca.Function("f", [S], [S - S]), [sv]),
      "times_zero": (ca.Function("f", [S], [S * 0]), [sv]),
      "full_reverse": (ca.Function("f", [A], [A[::-1, ::-1]]), [av]),
      "nz_duplicates": (ca.Function("f", [A], [A.nz[[1, 1, 2]]]), [av]),
      "single_element": (ca.Function("f", [A], [A[1, 2]]), [av]),
    }
    for name, (f, args) in cases.items():
      fn = "test_edge_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        g = ca.GraphBuilder(fn).create("imported_edge_%s" % name, {"symbolic": True})
        for o in range(f.n_out()):
          self.assertTrue(f.sparsity_out(o) == g.sparsity_out(o), "%s: out sparsity" % name)
        self.checkarray(ca.DM(f(*args)), ca.DM(g(*args)), name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def test_diagcat_blockdiag(self):
    # GENUINE block-diagonal assembly (diagcat outputs and intermediate OP_DIAGCAT) is done the
    # PSEUDO-DENSE way: each DENSE block is Pad-ded to the output shape at its (row,col) offset and
    # the padded blocks are Sum-med -- NO flatten/gather. A triangular block is just a dense block
    # with a structural zero, padded into place; the output sparsity is restored by the overlay.
    # PURE horzcat / vertcat outputs are not block-diagonal: they are a single ONNX Concat (no Pad).
    self.message("ONNX block assembly: Concat for pure horz/vertcat, Pad+Sum for diagcat")
    S = ca.MX.sym("S", ca.Sparsity.lower(3))
    sv = ca.DM(ca.Sparsity.lower(3), range(1, 7))
    a = ca.MX.sym("a", 3, 3); av = ca.DM(numpy.arange(9).reshape(3, 3).T)
    p, q = ca.MX.sym("p", 5, 2), ca.MX.sym("q", 5, 3)
    pv = ca.DM(numpy.arange(10).reshape(2, 5).T); qv = ca.DM(numpy.arange(15).reshape(3, 5).T)
    x, y = ca.MX.sym("x"), ca.MX.sym("y")
    cases = {
      "diag_sparse_block": (ca.Function("f", [S], [ca.diagcat(S, ca.MX.ones(2, 2))]), [sv]),
      "diag_dense_block": (ca.Function("f", [a], [ca.diagcat(a, ca.MX.ones(2, 2))]), [av]),
      "diag_three_blocks": (ca.Function("f", [S], [ca.diagcat(S, ca.MX.ones(2, 2), 7 * ca.MX.eye(2))]), [sv]),
      "horzcat": (ca.Function("f", [p, q], [ca.horzcat(p, q)]), [pv, qv]),
      "vertcat_scalars": (ca.Function("f", [x, y], [ca.vertcat(x, y)]), [3.0, 2.0]),
      "horzcat_scalars": (ca.Function("f", [x, y], [ca.horzcat(x, y)]), [3.0, 2.0]),
      "intermediate_diagcat": (ca.Function("f", [S], [ca.diagcat(S, S) * 2]), [sv]),
      "intermediate_then_slice": (ca.Function("f", [S], [(ca.diagcat(S, S) + 1)[4, :]]), [sv]),
      "intermediate_matmul": (ca.Function("f", [S], [ca.mtimes(ca.diagcat(S, ca.MX.ones(2, 2)),
                                                              ca.MX.ones(5, 1))]), [sv]),
    }
    # Genuine block-diagonal cases -> Pad+Sum. Pure horzcat/vertcat outputs -> a single Concat.
    padblock = {"diag_sparse_block", "diag_dense_block", "diag_three_blocks", "intermediate_diagcat"}
    concatblock = {"horzcat", "vertcat_scalars", "horzcat_scalars"}
    for name, (f, args) in cases.items():
      fn = "test_dc_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        if have_onnx and (name in padblock or name in concatblock):
          ops = [n.op_type for n in onnx.load(fn).graph.node]
          if name in padblock:
            self.assertTrue("Pad" in ops, "%s: expected Pad (block placement), got %s" % (name, ops))
          else:
            self.assertTrue("Concat" in ops, "%s: expected Concat (pure cat), got %s" % (name, ops))
            self.assertTrue("Pad" not in ops, "%s: pure cat should not Pad: %s" % (name, ops))
          # PSEUDO-DENSE: no flatten -- no Gather/ScatterElements/Reshape anywhere
          for bad in ("Gather", "ScatterElements", "Reshape"):
            self.assertTrue(bad not in ops, "%s: should not flatten (%s): %s" % (name, bad, ops))
        g = ca.GraphBuilder(fn).create("imported_dc_%s" % name, {"symbolic": True})
        for o in range(f.n_out()):
          self.assertTrue(f.sparsity_out(o) == g.sparsity_out(o), "%s: out sparsity" % name)
        self.checkarray(ca.DM(f(*args)), ca.DM(g(*args)), name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

  def _roundtrip_numeric(self, name, f, inputs):
    fn = "test_%s.onnx" % name
    try:
      ca.GraphBuilder(f).export_onnx(fn)
      g = ca.GraphBuilder(fn).create("imported_%s" % name, {"symbolic": True})
      for vals in inputs:
        vals = vals if isinstance(vals, tuple) else (vals,)
        self.checkarray(ca.DM(f(*vals)), ca.DM(g(*vals)), name, digits=12)
    finally:
      if os.path.exists(fn):
        os.remove(fn)

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

  def test_sparse_split(self):
    # horzsplit/vertsplit of a SPARSE value. The exporter must emit Split sizes in COLUMN/ROW units
    # (read off the segment dimensions), NOT nonzero-offset/stride (which only works for a dense
    # input). Each segment is a sub-block of the sparse seed, so its PATTERN round-trips natively.
    self.message("ONNX horzsplit/vertsplit of a sparse value (column/row split sizes)")
    # intricate 4x4 pattern via a sparse constant * scalar (a sparse seed), distinct values
    r = ca.DM.triplet([0, 1, 2, 2, 3, 1], [0, 0, 1, 3, 2, 3], [5, 7, 9, 11, 13, 15], 4, 4)
    s = ca.MX.sym("s")
    cases = {
      "hsplit2": ca.horzsplit(r * s, [0, 2, 4])[1],
      "hsplit_uneven": ca.horzsplit(r * s, [0, 1, 4])[1],
      "vsplit2": ca.vertsplit(r * s, [0, 2, 4])[0],
      "vsplit_uneven": ca.vertsplit(r * s, [0, 3, 4])[1],
    }
    for name, expr in cases.items():
      self.assertTrue(not expr.sparsity().is_dense(), "%s should be sparse" % name)
      f = ca.Function("f", [s], [expr])
      fn = "test_ssplit_%s.onnx" % name
      try:
        ca.GraphBuilder(f).export_onnx(fn)
        g = ca.GraphBuilder(fn).create("imported_ssplit_%s" % name, {"symbolic": True})
        self.assertTrue(f.sparsity_out(0) == g.sparsity_out(0), "%s: sparsity not preserved" % name)
        for v in [2.0, -1.5]:
          self.checkarray(f(v), g(v), "ssplit_%s" % name, digits=12)
      finally:
        if os.path.exists(fn):
          os.remove(fn)

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
    # Function.__call__ is typed -> MX for the common single-output case; base has 2 outputs so it returns a list
    self.check("map_multi", ca.Function("f", [X, Y], base.map(3)(X, Y)),  # pyright: ignore[reportArgumentType,reportCallIssue]
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

  def test_kron(self):
    # Kronecker product kron(A,B). Export lowers (stored space, no extra transposes since kron
    # commutes with transpose) to Reshape(A->[ca,1,ra,1]) * Reshape(B->[1,cb,1,rb]) broadcast,
    # then Reshape back; the node group is name-tagged "kron<k>_*" so import reconstructs
    # kron(A_imp,B_imp) directly (the 4-D intermediates cannot be 2-D MX). Sparsity round-trips
    # natively (kron propagates patterns).
    self.message("ONNX kron (Kronecker product) round-trip")

    # Dense cases -> full ORT validation via check()
    A = ca.MX.sym("A", 2, 3)
    B = ca.MX.sym("B", 2, 2)
    self.check("kron_dense_2x3_2x2", ca.Function("f", [A, B], [ca.kron(A, B)]),
               [(ca.DM([[1, 2, 3], [4, 5, 6]]), ca.DM([[1, 2], [3, 4]]))])

    u = ca.MX.sym("u", 3, 1)
    v = ca.MX.sym("v", 1, 2)
    self.check("kron_dense_vec", ca.Function("f", [u, v], [ca.kron(u, v)]),
               [(ca.DM([[1], [2], [3]]), ca.DM([[10, 20]]))])

    C = ca.MX.sym("C", 1, 1)
    D = ca.MX.sym("D", 3, 2)
    self.check("kron_dense_scalar", ca.Function("f", [C, D], [ca.kron(C, D)]),
               [(ca.DM(5), ca.DM([[1, 2], [3, 4], [5, 6]]))])

    # Sparse-seed case: the sparse function input plants an all-zero sparse_initializer; import
    # recovers the pattern and kron propagates it. (ORT shadows a sparse-initializer input out of
    # get_inputs(), so the standard positional-feed onnxruntime_check can't drive it; we validate
    # ORT explicitly below by feeding all graph inputs by name.)
    S = ca.MX.sym("S", ca.Sparsity.lower(3))
    T = ca.MX.sym("T", 2, 2)
    fs = ca.Function("f", [S, T], [ca.kron(S, T)])
    Sv = ca.DM(ca.Sparsity.lower(3), [1, 2, 3, 4, 5, 6])
    Tv = ca.DM([[1, 0], [0, 2]])
    fn = "test_kron_sparse.onnx"
    try:
      ca.GraphBuilder(fs).export_onnx(fn)
      g = ca.GraphBuilder(fn).create("imported_kron_sparse", {"symbolic": True})
      # Sparsity preserved (not densified) and numerically identical
      self.assertTrue(fs.sparsity_out(0) == g.sparsity_out(0), "kron sparsity not preserved")
      self.assertEqual(g.nnz_out(0), fs(Sv, Tv).nnz())
      self.compare(fs(Sv, Tv), g(Sv, Tv), "kron_sparse")
      if has_onnxruntime:
        session = ort.InferenceSession(fn)
        # Feed every graph input by name (sparse-initializer-shadowed inputs are omitted from
        # get_inputs() but remain overridable by name); un-transpose ORT's output.
        feed = {n: numpy.array(ca.densify(ca.DM(v))).T.astype(numpy.float64)
                for n, v in zip(fs.name_in(), (Sv, Tv))}
        out = session.run(None, feed)[0]
        self.checkarray(fs(Sv, Tv), ca.DM(numpy.array(out).T), "kron_sparse_ort", digits=12)
    finally:
      if os.path.exists(fn):
        os.remove(fn)

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
