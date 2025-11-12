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
from casadi import *

import casadi as ca

import os
import unittest
from helpers import *
import numpy as np

class Onnxtests(casadiTestCase):

  def validate_with_onnxruntime(self, onnx_file, casadi_function, test_inputs):
    """Validate CasADi function output against ONNX Runtime"""
    try:
      import onnxruntime as ort
    except ImportError:
      print(f"onnxruntime not available - skipping numerical validation for {onnx_file}")
      return

    # Load ONNX model with ONNX Runtime
    session = ort.InferenceSession(onnx_file)

    # Get input/output names
    input_names = [inp.name for inp in session.get_inputs()]

    # Prepare inputs for ONNX Runtime
    if not isinstance(test_inputs, list):
      test_inputs = [test_inputs]

    onnx_inputs = {}
    for i, inp_name in enumerate(input_names):
      casadi_inp = test_inputs[i] if i < len(test_inputs) else DM.zeros(1)
      onnx_inputs[inp_name] = np.array(casadi_inp).astype(np.float64)

    # Run ONNX Runtime inference
    onnx_outputs = session.run(None, onnx_inputs)

    # Run CasADi function
    casadi_outputs = casadi_function(*test_inputs)
    if not isinstance(casadi_outputs, list):
      casadi_outputs = [casadi_outputs]

    # Compare outputs
    for i, (onnx_out, casadi_out) in enumerate(zip(onnx_outputs, casadi_outputs)):
      self.checkarray(DM(onnx_out), casadi_out, digits=10)

  def test_export_numerical_validation(self):
    """Test ONNX export with numerical validation against ONNX Runtime"""
    if not has_translator("onnx"):
      self.skipTest("ONNX translator not available")

    # Create a function with operations we support: sin, mul, add
    x = MX.sym("x")
    y = MX.sym("y")
    z = sin(x) * y + x
    f = Function("test_numerical", [x, y], [z], ["x", "y"], ["z"])

    # Export to ONNX
    onnx_file = "test_onnx_numerical.onnx"
    t = translator("onnx")
    t.load(f)
    t.save(onnx_file)

    # Check file was created
    self.assertTrue(os.path.exists(onnx_file))

    # Validate against ONNX Runtime if available
    test_inputs = [DM(0.5), DM(1.2)]
    self.validate_with_onnxruntime(onnx_file, f, test_inputs)

    # Cleanup
    if os.path.exists(onnx_file):
      os.remove(onnx_file)

  def test_import_roundtrip(self):
    """Test ONNX import via roundtrip (export then import)"""
    if not has_translator("onnx"):
      self.skipTest("ONNX translator not available")

    # Create a simple function using supported operations: sin, mul, add
    x = MX.sym("x")
    y = MX.sym("y")
    z = sin(x) * y + x
    f_original = Function("test_roundtrip", [x, y], [z], ["x", "y"], ["z"])

    # Export to ONNX
    onnx_file = "test_roundtrip.onnx"
    t_export = translator("onnx", {"verbose": False})
    t_export.load(f_original)
    t_export.save(onnx_file)

    # Verify file was created
    self.assertTrue(os.path.exists(onnx_file), "ONNX file should be created")

    # Import from ONNX
    t_import = translator("onnx", {"verbose": False})
    t_import.load(onnx_file)
    f_imported = t_import.create("imported_function")

    # Verify function properties
    self.assertEqual(f_imported.n_in(), 2, "Should have 2 inputs")
    self.assertEqual(f_imported.n_out(), 1, "Should have 1 output")
    self.assertEqual(f_imported.name_in(0), "x", "First input should be 'x'")
    self.assertEqual(f_imported.name_in(1), "y", "Second input should be 'y'")
    self.assertEqual(f_imported.name_out(0), "z", "Output should be 'z'")

    # Numerical validation - test several input combinations
    test_cases = [
      [DM(0.5), DM(1.2)],
      [DM(0.0), DM(1.0)],
      [DM(1.5), DM(2.3)],
      [DM(-0.5), DM(0.8)]
    ]

    for test_inputs in test_cases:
      result_original = f_original(*test_inputs)
      result_imported = f_imported(*test_inputs)

      # Compare outputs (should match exactly within floating point precision)
      self.checkarray(result_original, result_imported, digits=10)

    # Cleanup
    if os.path.exists(onnx_file):
      os.remove(onnx_file)

  def test_import_all_operations(self):
    """Test ONNX import with all 14 supported operations"""
    if not has_translator("onnx"):
      self.skipTest("ONNX translator not available")

    # Create a function using all supported operations
    x = MX.sym("x")
    y = MX.sym("y")

    # Test all operations:
    # Binary: Add, Sub, Mul, Div
    # Unary: Sin, Cos, Tan, Exp, Log, Sqrt, Neg, Tanh
    z = (sin(x) + cos(x) - tan(y) * exp(y)) / (log(x + 2) + sqrt(y + 1)) + tanh(x) - (-y)

    f_original = Function("test_all_ops", [x, y], [z], ["x", "y"], ["z"])

    # Export to ONNX
    onnx_file = "test_all_ops.onnx"
    t_export = translator("onnx", {"verbose": False})
    t_export.load(f_original)
    t_export.save(onnx_file)

    # Import from ONNX
    t_import = translator("onnx", {"verbose": False})
    t_import.load(onnx_file)
    f_imported = t_import.create("imported_all_ops")

    # Numerical validation with multiple test cases
    test_cases = [
      [DM(1.0), DM(0.5)],
      [DM(2.0), DM(1.0)],
      [DM(0.5), DM(0.2)],
      [DM(1.5), DM(1.5)]
    ]

    for test_inputs in test_cases:
      result_original = f_original(*test_inputs)
      result_imported = f_imported(*test_inputs)

      # Compare outputs
      self.checkarray(result_original, result_imported, digits=10)

    # Cleanup
    if os.path.exists(onnx_file):
      os.remove(onnx_file)

if __name__ == '__main__':
    unittest.main()
