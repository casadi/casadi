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

import unittest
import os
from casadi import *
from casadi.tools import *

# Import base test class
from helpers import casadiTestCase

# ============================================================================
# GLOBAL AVAILABILITY FLAGS
# ============================================================================

# Check ONNX Runtime availability once at module load
try:
    import onnxruntime as ort
    import numpy as np
    HAS_ONNXRUNTIME = True
except ImportError:
    HAS_ONNXRUNTIME = False

# Tests to skip for ONNX Runtime validation (test_name -> reason)
SKIP_ONNXRUNTIME = {
    # Float64 not supported in ONNX Runtime for these operations
    'unary_cos': "No float64 support",
    'unary_tan': "No float64 support",
    'unary_asin': "No float64 support",
    'unary_acos': "No float64 support",
    'unary_atan': "No float64 support",
    'unary_sinh': "No float64 support",
    'unary_cosh': "No float64 support",
    'unary_asinh': "No float64 support",
    'unary_acosh': "No float64 support",
    'unary_atanh': "No float64 support",
    'unary_erf': "No float64 support",
    # Empty function has shape mismatch
    'empty_function': "Shape mismatch in ONNX Runtime",
    # Boolean/comparison ops: ONNX requires boolean tensor types
    'logic_not': "ONNX Not requires boolean input",
    'logic_or': "ONNX Or requires boolean inputs",
    'logic_and': "ONNX And requires boolean inputs",
    'ne': "ONNX Equal/Not output boolean, graph declares double",
    'eq': "ONNX Equal outputs boolean, graph declares double",
    'le': "ONNX LessOrEqual outputs boolean, graph declares double",
    'lt': "ONNX Less outputs boolean, graph declares double",
    # Transpose+MatMul fusion: ONNX Runtime fuses to FusedMatMul which lacks float64 support
    'bilin': "FusedMatMul optimization lacks float64 support",
    'rank1': "FusedMatMul optimization lacks float64 support",
}

# ============================================================================
# TEST DATA
# ============================================================================

# Unary operations: (name, expr_lambda, test_values)
UNARY_OPS = [
    ("sin", lambda x: sin(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("cos", lambda x: cos(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("tan", lambda x: tan(x), [DM(0.0), DM(0.5), DM(1.0)]),
    ("exp", lambda x: exp(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)]),
    ("log", lambda x: log(x), [DM(0.1), DM(1.0), DM(2.0), DM(10.0)]),
    ("sqrt", lambda x: sqrt(x), [DM(0.0), DM(1.0), DM(4.0), DM(9.0)]),
    ("asin", lambda x: asin(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("acos", lambda x: acos(x), [DM(0.0), DM(0.5), DM(1.0)]),
    ("atan", lambda x: atan(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)]),
    ("sinh", lambda x: sinh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("cosh", lambda x: cosh(x), [DM(0.0), DM(0.5), DM(1.0)]),
    ("tanh", lambda x: tanh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)]),
    ("asinh", lambda x: asinh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("acosh", lambda x: acosh(x), [DM(1.0), DM(1.5), DM(2.0)]),
    ("atanh", lambda x: atanh(x), [DM(0.0), DM(0.5), DM(-0.5)]),
    ("ceil", lambda x: ceil(x), [DM(1.2), DM(1.5), DM(1.9), DM(-1.2)]),
    ("floor", lambda x: floor(x), [DM(1.2), DM(1.5), DM(1.9), DM(-1.2)]),
    ("fabs", lambda x: fabs(x), [DM(1.0), DM(-1.0), DM(0.0), DM(-5.5)]),
    ("sign", lambda x: sign(x), [DM(1.0), DM(-1.0), DM(0.0), DM(5.5)]),
    ("neg", lambda x: -x, [DM(1.0), DM(-1.0), DM(0.0), DM(5.5)]),
    ("erf", lambda x: erf(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
]

# Binary operations: (name, expr_lambda, test_values)
BINARY_OPS = [
    ("add", lambda x, y: x + y, [(DM(2.0), DM(3.0)), (DM(-1.0), DM(4.0)), (DM(0.0), DM(0.0))]),
    ("sub", lambda x, y: x - y, [(DM(5.5), DM(2.2)), (DM(10.0), DM(3.0)), (DM(0.0), DM(0.0))]),
    ("mul", lambda x, y: x * y, [(DM(2.0), DM(3.0)), (DM(-2.0), DM(4.0)), (DM(0.0), DM(5.0))]),
    ("div", lambda x, y: x / y, [(DM(10.0), DM(2.0)), (DM(7.0), DM(3.0)), (DM(1.0), DM(1.0))]),
    ("pow", lambda x, y: x ** y, [(DM(2.0), DM(3.0)), (DM(4.0), DM(0.5)), (DM(10.0), DM(2.0))]),
]

# Matrix operations: (name, expr_lambda, shapes, test_values)
MATRIX_OPS = [
    ("matmul_2x3_3x2", lambda A, B: mtimes(A, B), [(2, 3), (3, 2)],
     [(DM([[1, 2, 3], [4, 5, 6]]), DM([[1, 2], [3, 4], [5, 6]]))]),
    ("matmul_3x1_1x3", lambda A, B: mtimes(A, B), [(3, 1), (1, 3)],
     [(DM([[1], [2], [3]]), DM([[4, 5, 6]]))]),
    ("transpose_2x3", lambda A: A.T, [(2, 3)],
     [(DM([[1, 2, 3], [4, 5, 6]]),)]),
]


# ============================================================================
# TEST CLASS
# ============================================================================

class Onnxtests(casadiTestCase):
    """ONNX translator tests"""

    def roundtrip_test(self, op_name, casadi_func, test_inputs):
        """
        Execute roundtrip test: CasADi -> ONNX -> CasADi -> validate

        Args:
            op_name: Name of operation (for file naming)
            casadi_func: CasADi Function to test
            test_inputs: List of input values to test (each can be single value or tuple)
        """
        onnx_file = f"test_{op_name}.onnx"

        try:
            # Export to ONNX
            t_export = translator("onnx", {"verbose": False})
            t_export.load(casadi_func)
            t_export.save(onnx_file)

            # Import from ONNX
            t_import = translator("onnx", {"verbose": False})
            t_import.load(onnx_file)
            f_imported = t_import.create(f"imported_{op_name}")

            # Validate numerical results
            for test_val in test_inputs:
                # Handle both single values and tuples
                if not isinstance(test_val, tuple):
                    test_val = (test_val,)

                result_original = casadi_func(*test_val)
                result_imported = f_imported(*test_val)

                # Handle different return types
                if isinstance(result_original, dict):
                    for key in result_original:
                        self.checkarray(result_original[key], result_imported[key], digits=10)
                elif isinstance(result_original, (list, tuple)):
                    for orig, imp in zip(result_original, result_imported):
                        self.checkarray(orig, imp, digits=10)
                else:
                    self.checkarray(result_original, result_imported, digits=10)

        finally:
            # Cleanup
            if os.path.exists(onnx_file):
                os.remove(onnx_file)

    def onnxruntime_test(self, op_name, casadi_func, test_inputs):
        """
        Test export and validate with ONNX Runtime

        Args:
            op_name: Name of operation
            casadi_func: CasADi Function to test
            test_inputs: Input values for testing
        """
        onnx_file = f"test_{op_name}_onnxruntime.onnx"

        try:
            # Export
            t = translator("onnx", {"verbose": False})
            t.load(casadi_func)
            t.save(onnx_file)

            # Load ONNX model and get input names from session
            session = ort.InferenceSession(onnx_file)
            onnx_input_names = [inp.name for inp in session.get_inputs()]

            # Validate all test inputs
            for test_val in test_inputs:
                # Handle both single values and tuples
                if not isinstance(test_val, tuple):
                    test_val = (test_val,)

                # Convert to numpy format using ONNX model's input names
                onnx_inputs = {}
                for i, inp in enumerate(test_val):
                    input_name = onnx_input_names[i] if i < len(onnx_input_names) else f"input_{i}"
                    onnx_inputs[input_name] = np.array(inp).astype(np.float64)

                # Run CasADi
                casadi_outputs = casadi_func(*test_val)
                if not isinstance(casadi_outputs, (list, tuple)):
                    casadi_outputs = [casadi_outputs]

                # Run ONNX Runtime
                onnx_outputs = session.run(None, onnx_inputs)

                # Compare
                for i, (onnx_out, casadi_out) in enumerate(zip(onnx_outputs, casadi_outputs)):
                    self.checkarray(DM(onnx_out), casadi_out, digits=10)

        finally:
            if os.path.exists(onnx_file):
                os.remove(onnx_file)

    @classmethod
    def add_test(cls, test_name, func, test_inputs, doc=None):
        """
        Add a test (roundtrip + onnxruntime if available) to the test class.

        Args:
            test_name: Base name for the test (e.g., "unary_sin")
            func: CasADi Function to test
            test_inputs: List of test input values
            doc: Optional docstring for the test
        """
        # Create roundtrip test
        def make_roundtrip_test(name, f, inputs):
            def test_method(self):
                self.roundtrip_test(name, f, inputs)
            return test_method

        roundtrip_test = make_roundtrip_test(test_name, func, test_inputs)
        roundtrip_test.__name__ = f"test_{test_name}_roundtrip"
        roundtrip_test.__doc__ = doc or f"Test {test_name} via roundtrip"
        setattr(cls, roundtrip_test.__name__, roundtrip_test)

        # Create onnxruntime test
        def make_onnxruntime_test(name, f, inputs):
            def test_method(self):
                self.onnxruntime_test(name, f, inputs)
            return test_method

        onnxruntime_test = make_onnxruntime_test(test_name, func, test_inputs)
        onnxruntime_test.__name__ = f"test_{test_name}_onnxruntime"
        onnxruntime_test.__doc__ = doc or f"Test {test_name} via onnxruntime"

        # Apply skip decorator based on SKIP_ONNXRUNTIME dict or availability
        if test_name in SKIP_ONNXRUNTIME:
            onnxruntime_test = unittest.skip(SKIP_ONNXRUNTIME[test_name])(onnxruntime_test)
        elif not HAS_ONNXRUNTIME:
            onnxruntime_test = unittest.skip("ONNX Runtime not available")(onnxruntime_test)

        setattr(cls, onnxruntime_test.__name__, onnxruntime_test)


# ============================================================================
# TEST REGISTRATION
# ============================================================================

# Shared symbols for test registration
_x = MX.sym("x")
_y = MX.sym("y")
_z = MX.sym("z")


def _register_op_tests():
    """Register data-driven tests for operations"""

    # Unary operations
    for op_name, expr, test_values in UNARY_OPS:
        Onnxtests.add_test(
            f"unary_{op_name}",
            Function(f"test_{op_name}", [_x], [expr(_x)]),
            test_values,
            doc=f"Test unary {op_name} operation"
        )

    # Binary operations
    for op_name, expr, test_values in BINARY_OPS:
        Onnxtests.add_test(
            f"binary_{op_name}",
            Function(f"test_{op_name}", [_x, _y], [expr(_x, _y)]),
            test_values,
            doc=f"Test binary {op_name} operation"
        )

    # Matrix operations
    for op_name, expr, shapes, test_values in MATRIX_OPS:
        if len(shapes) == 1:
            A = MX.sym("A", shapes[0][0], shapes[0][1])
            f = Function(f"test_{op_name}", [A], [expr(A)])
        else:
            A = MX.sym("A", shapes[0][0], shapes[0][1])
            B = MX.sym("B", shapes[1][0], shapes[1][1])
            f = Function(f"test_{op_name}", [A, B], [expr(A, B)])

        Onnxtests.add_test(
            f"matrix_{op_name}",
            f,
            test_values,
            doc=f"Test matrix {op_name} operation"
        )


def _register_complex_tests():
    """Register tests for complex expressions, control flow, and function hierarchies"""

    # Different input shapes for sin function
    shapes_to_test = [
        ((1, 1), DM([[2.0]])),
        ((3, 1), DM([[1.0], [2.0], [3.0]])),
        ((1, 3), DM([[1.0, 2.0, 3.0]])),
        ((2, 2), DM([[1.0, 2.0], [3.0, 4.0]])),
    ]

    for shape, test_val in shapes_to_test:
        x = MX.sym("x", shape[0], shape[1])
        f = Function(f"sin_{shape[0]}x{shape[1]}", [x], [sin(x)])
        Onnxtests.add_test(f"sin_{shape[0]}x{shape[1]}", f, [test_val])

    # Empty function (no inputs)
    Onnxtests.add_test(
        "empty_function",
        Function("empty", [], [MX(42.0)]),
        [()],
        doc="Test function with no inputs"
    )

    # Complex expression (using sin instead of cos for ONNX Runtime float64 compatibility)
    Onnxtests.add_test(
        "complex_expression",
        Function("complex_expr", [_x, _y], [(sin(_x) + sin(_y)) * exp(_x - _y)]),
        [(DM(0.5), DM(1.0)), (DM(1.0), DM(0.5)), (DM(0.0), DM(0.0))],
        doc="Test complex expression combining multiple operations"
    )

    # Multiple outputs
    Onnxtests.add_test(
        "multiple_outputs",
        Function("multi_out", [_x, _y], [_x + _y, _x - _y, _x * _y]),
        [(DM(3.0), DM(2.0)), (DM(5.0), DM(1.0))],
        doc="Test function with multiple outputs"
    )

    # Function hierarchy: simple
    f_inner = Function("inner", [_x], [_x + _x])
    Onnxtests.add_test(
        "function_hierarchy_simple",
        Function("outer", [_y], [f_inner(_y) + 1]).wrap(),
        [DM(5.0)],
        doc="Test nested function call: outer calls inner"
    )

    # Function hierarchy: multiple calls
    f_double = Function("double", [_x], [_x + _x])
    Onnxtests.add_test(
        "function_hierarchy_multiple_calls",
        Function("outer_multi", [_y], [f_double(_y) + f_double(_y + 1)]).wrap(),
        [DM(3.0)],
        doc="Test outer function calling inner function twice"
    )

    # Function hierarchy: deep (A calls B calls C)
    f_c = Function("func_c", [_x], [sin(_x)])
    f_b = Function("func_b", [_y], [f_c(_y) + 1])
    Onnxtests.add_test(
        "function_hierarchy_deep",
        Function("func_a", [_z], [f_b(_z) * 2]).wrap(),
        [DM(0.5)],
        doc="Test deep function hierarchy: A calls B calls C"
    )

    # Vertcat input
    x = MX.sym('x')
    y = MX.sym('y')
    Onnxtests.add_test(
        "vertcat_input",
        Function("vertcat_test", [vertcat(x, y)], [x + y]),
        [DM([3.0, 2.0])],  # Input is [x; y] = [3; 2], output is x+y = 5
        doc="Test vertcat input"
    )

    # Vertcat output
    x = MX.sym('x')
    y = MX.sym('y')
    Onnxtests.add_test(
        "vertcat_output",
        Function("vertcat_out", [x, y], [vertcat(x, y)]),
        [(DM(3.0), DM(2.0))],  # Inputs x=3, y=2, output is [3; 2]
        doc="Test vertcat output"
    )


def _register_mx_ops_tests():
    """Register tests for MX operations (repmat, indexing, etc.)"""

    # ======= Repmat tests =======
    # repmat horizontal (1, n)
    x_2x2 = MX.sym("x", 2, 2)
    Onnxtests.add_test(
        "repmat_1x3",
        Function("repmat_1x3", [x_2x2], [repmat(x_2x2, 1, 3)]),
        [DM([[1, 2], [3, 4]])],
        doc="Test repmat with 1 row, 3 column repetitions"
    )

    # repmat with n=2
    x_3x1 = MX.sym("x", 3, 1)
    Onnxtests.add_test(
        "repmat_1x2",
        Function("repmat_1x2", [x_3x1], [repmat(x_3x1, 1, 2)]),
        [DM([[1], [2], [3]])],
        doc="Test repmat with 1 row, 2 column repetitions"
    )

    # ======= Indexing tests =======
    # Single element indexing from vector - x[0]
    x_vec = MX.sym("x", 5, 1)
    Onnxtests.add_test(
        "index_single_element",
        Function("idx_single", [x_vec], [x_vec[0]]),
        [DM([1, 2, 3, 4, 5])],
        doc="Test single element indexing x[0]"
    )

    # Multiple element indexing - x[[0,2,4]]
    x_vec5 = MX.sym("x", 5, 1)
    Onnxtests.add_test(
        "index_multiple_elements",
        Function("idx_multi", [x_vec5], [x_vec5[[0, 2, 4]]]),
        [DM([10, 20, 30, 40, 50])],
        doc="Test multiple element indexing x[[0,2,4]]"
    )

    # Slice indexing - x[1:4]
    x_vec6 = MX.sym("x", 6, 1)
    Onnxtests.add_test(
        "index_slice",
        Function("idx_slice", [x_vec6], [x_vec6[1:4]]),
        [DM([1, 2, 3, 4, 5, 6])],
        doc="Test slice indexing x[1:4]"
    )

    # ======= Dot product test =======
    x_dot = MX.sym("x", 3, 1)
    y_dot = MX.sym("y", 3, 1)
    Onnxtests.add_test(
        "dot_product",
        Function("dot_prod", [x_dot, y_dot], [dot(x_dot, y_dot)]),
        [(DM([1, 2, 3]), DM([4, 5, 6]))],
        doc="Test dot product operation"
    )

    # ======= Repsum test =======
    x_repsum = MX.sym("x", 2, 6)
    Onnxtests.add_test(
        "repsum_1x3",
        Function("repsum_1x3", [x_repsum], [repsum(x_repsum, 1, 3)]),
        [DM([[1, 2, 10, 20, 100, 200], [3, 4, 30, 40, 300, 400]])],
        doc="Test repsum with 1 row, 3 column sums"
    )

    # ======= Blockcat test =======
    a_bc = MX.sym("a", 2, 2)
    b_bc = MX.sym("b", 2, 2)
    c_bc = MX.sym("c", 2, 2)
    d_bc = MX.sym("d", 2, 2)
    Onnxtests.add_test(
        "blockcat_2x2",
        Function("blockcat_2x2", [a_bc, b_bc, c_bc, d_bc], [blockcat([[a_bc, b_bc], [c_bc, d_bc]])]),
        [(DM([[1, 2], [3, 4]]), DM([[5, 6], [7, 8]]), DM([[9, 10], [11, 12]]), DM([[13, 14], [15, 16]]))],
        doc="Test blockcat 2x2 block matrix"
    )


def _register_tensor_op_tests():
    """Register tests for tensor operations (split, concat, slice, reduce)"""

    # ======= HorzSplit tests =======
    # Split a 2x6 matrix into three 2x2 parts, use them
    x_hs = MX.sym("x", 2, 6)
    parts = horzsplit(x_hs, [0, 2, 4, 6])
    Onnxtests.add_test(
        "horzsplit_equal",
        Function("horzsplit_eq", [x_hs], [parts[0] + parts[1] + parts[2]]),
        [DM([[1, 2, 10, 20, 100, 200], [3, 4, 30, 40, 300, 400]])],
        doc="Test horzsplit with equal-size splits"
    )

    # Unequal split: 3x5 -> 3x2 and 3x3
    x_hs2 = MX.sym("x", 3, 5)
    parts2 = horzsplit(x_hs2, [0, 2, 5])
    Onnxtests.add_test(
        "horzsplit_unequal",
        Function("horzsplit_uneq", [x_hs2], [parts2[0], parts2[1]]),
        [DM([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]])],
        doc="Test horzsplit with unequal splits"
    )

    # ======= VertSplit tests =======
    # Split a 6x2 matrix into three 2x2 parts
    x_vs = MX.sym("x", 6, 2)
    vparts = vertsplit(x_vs, [0, 2, 4, 6])
    Onnxtests.add_test(
        "vertsplit_equal",
        Function("vertsplit_eq", [x_vs], [vparts[0] + vparts[1] + vparts[2]]),
        [DM([[1, 2], [3, 4], [10, 20], [30, 40], [100, 200], [300, 400]])],
        doc="Test vertsplit with equal-size splits"
    )

    # Unequal split: 5x2 -> 2x2 and 3x2
    x_vs2 = MX.sym("x", 5, 2)
    vparts2 = vertsplit(x_vs2, [0, 2, 5])
    Onnxtests.add_test(
        "vertsplit_unequal",
        Function("vertsplit_uneq", [x_vs2], [vparts2[0], vparts2[1]]),
        [DM([[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]])],
        doc="Test vertsplit with unequal splits"
    )

    # ======= Concat tests =======
    # Horizontal concat: two 2x2 -> 2x4
    a_hc = MX.sym("a", 2, 2)
    b_hc = MX.sym("b", 2, 3)
    Onnxtests.add_test(
        "horzcat_2x2_2x3",
        Function("horzcat_test", [a_hc, b_hc], [horzcat(a_hc, b_hc)]),
        [(DM([[1, 2], [3, 4]]), DM([[5, 6, 7], [8, 9, 10]]))],
        doc="Test horizontal concatenation"
    )

    # Vertical concat: 2x3 and 3x3 -> 5x3
    a_vc = MX.sym("a", 2, 3)
    b_vc = MX.sym("b", 3, 3)
    Onnxtests.add_test(
        "vertcat_2x3_3x3",
        Function("vertcat_test", [a_vc, b_vc], [vertcat(a_vc, b_vc)]),
        [(DM([[1, 2, 3], [4, 5, 6]]), DM([[7, 8, 9], [10, 11, 12], [13, 14, 15]]))],
        doc="Test vertical concatenation"
    )

    # ======= Slice tests =======
    # 2D slice: extract submatrix
    x_sl = MX.sym("x", 4, 4)
    Onnxtests.add_test(
        "slice_2d",
        Function("slice_2d", [x_sl], [x_sl[1:3, 0:2]]),
        [DM([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])],
        doc="Test 2D slice x[1:3, 0:2]"
    )

    # Row slice
    x_rs = MX.sym("x", 4, 3)
    Onnxtests.add_test(
        "slice_rows",
        Function("slice_rows", [x_rs], [x_rs[1:3, :]]),
        [DM([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])],
        doc="Test row slice x[1:3, :]"
    )

    # ======= ReduceMin/ReduceMax tests =======
    x_rd = MX.sym("x", 3, 2)
    Onnxtests.add_test(
        "reduce_min",
        Function("reduce_min", [x_rd], [mmin(x_rd)]),
        [DM([[5, 2], [1, 8], [3, 4]])],
        doc="Test mmin (ReduceMin) operation"
    )

    Onnxtests.add_test(
        "reduce_max",
        Function("reduce_max", [x_rd], [mmax(x_rd)]),
        [DM([[5, 2], [1, 8], [3, 4]])],
        doc="Test mmax (ReduceMax) operation"
    )

    # ======= Reciprocal test =======
    x_inv = MX.sym("x", 2, 2)
    Onnxtests.add_test(
        "reciprocal",
        Function("reciprocal", [x_inv], [1.0 / x_inv]),
        [DM([[1, 2], [4, 5]])],
        doc="Test reciprocal (1/x) operation"
    )

    # ======= Norm tests =======
    x_n = MX.sym("x", 3, 1)
    Onnxtests.add_test(
        "norm_1",
        Function("norm_1", [x_n], [norm_1(x_n)]),
        [DM([1, -2, 3])],
        doc="Test L1 norm (ReduceL1) operation"
    )

    Onnxtests.add_test(
        "norm_2",
        Function("norm_2", [x_n], [norm_2(x_n)]),
        [DM([3, 4, 0])],
        doc="Test L2 norm (ReduceL2) operation"
    )


def _register_additional_tests():
    """Register tests for additional operations (sq, twice, logic, normf, ne, norminf)"""

    # ======= SQ test (x^2) =======
    x_sq = MX.sym("x", 2, 2)
    Onnxtests.add_test(
        "sq",
        Function("sq_test", [x_sq], [x_sq**2]),
        [DM([[1, 2], [3, 4]])],
        doc="Test square operation (x^2)"
    )

    # ======= TWICE test (2*x) =======
    x_tw = MX.sym("x", 2, 2)
    Onnxtests.add_test(
        "twice",
        Function("twice_test", [x_tw], [x_tw + x_tw]),
        [DM([[1, 2], [3, 4]])],
        doc="Test twice operation (2*x)"
    )

    # ======= NORMF test (Frobenius norm on matrix) =======
    x_nf = MX.sym("x", 3, 2)
    Onnxtests.add_test(
        "norm_fro",
        Function("norm_fro_test", [x_nf], [norm_fro(x_nf)]),
        [DM([[1, 2], [3, 4], [5, 6]])],
        doc="Test Frobenius norm (ReduceL2) on matrix"
    )

    # ======= NORMINF test (infinity norm) =======
    x_ni = MX.sym("x", 3, 1)
    Onnxtests.add_test(
        "norm_inf",
        Function("norm_inf_test", [x_ni], [norm_inf(x_ni)]),
        [DM([1, -5, 3])],
        doc="Test infinity norm (max(abs(x))) operation"
    )

    # ======= NE test (not equal) =======
    x_ne = MX.sym("x")
    y_ne = MX.sym("y")
    Onnxtests.add_test(
        "ne",
        Function("ne_test", [x_ne, y_ne], [x_ne != y_ne]),
        [(DM(1.0), DM(1.0)), (DM(1.0), DM(2.0))],
        doc="Test not-equal comparison"
    )

    # ======= Logic NOT test =======
    x_not = MX.sym("x")
    Onnxtests.add_test(
        "logic_not",
        Function("not_test", [x_not], [logic_not(x_not)]),
        [DM(0.0), DM(1.0)],
        doc="Test logical NOT operation"
    )

    # ======= Logic OR test =======
    x_or = MX.sym("x")
    y_or = MX.sym("y")
    Onnxtests.add_test(
        "logic_or",
        Function("or_test", [x_or, y_or], [logic_or(x_or, y_or)]),
        [(DM(0.0), DM(0.0)), (DM(0.0), DM(1.0)), (DM(1.0), DM(1.0))],
        doc="Test logical OR operation"
    )

    # ======= Logic AND test =======
    x_and = MX.sym("x")
    y_and = MX.sym("y")
    Onnxtests.add_test(
        "logic_and",
        Function("and_test", [x_and, y_and], [logic_and(x_and, y_and)]),
        [(DM(0.0), DM(0.0)), (DM(0.0), DM(1.0)), (DM(1.0), DM(1.0))],
        doc="Test logical AND operation"
    )

    # ======= FMOD test (modulo) =======
    x_fm = MX.sym("x", 2, 2)
    y_fm = MX.sym("y", 2, 2)
    Onnxtests.add_test(
        "fmod",
        Function("fmod_test", [x_fm, y_fm], [fmod(x_fm, y_fm)]),
        [(DM([[5, 7], [10, 3]]), DM([[3, 4], [3, 2]]))],
        doc="Test fmod (modulo) operation"
    )

    # ======= COPYSIGN test =======
    x_cs = MX.sym("x", 2, 2)
    y_cs = MX.sym("y", 2, 2)
    Onnxtests.add_test(
        "copysign",
        Function("copysign_test", [x_cs, y_cs], [copysign(x_cs, y_cs)]),
        [(DM([[5, -3], [-7, 2]]), DM([[-1, 1], [1, -1]]))],
        doc="Test copysign operation"
    )

    # ======= EQ test (equal) =======
    x_eq = MX.sym("x")
    y_eq = MX.sym("y")
    Onnxtests.add_test(
        "eq",
        Function("eq_test", [x_eq, y_eq], [x_eq == y_eq]),
        [(DM(1.0), DM(1.0)), (DM(1.0), DM(2.0))],
        doc="Test equality comparison"
    )

    # ======= LE test (less or equal) =======
    x_le = MX.sym("x")
    y_le = MX.sym("y")
    Onnxtests.add_test(
        "le",
        Function("le_test", [x_le, y_le], [x_le <= y_le]),
        [(DM(1.0), DM(2.0)), (DM(2.0), DM(2.0)), (DM(3.0), DM(2.0))],
        doc="Test less-or-equal comparison"
    )

    # ======= LT test (less than) =======
    x_lt = MX.sym("x")
    y_lt = MX.sym("y")
    Onnxtests.add_test(
        "lt",
        Function("lt_test", [x_lt, y_lt], [x_lt < y_lt]),
        [(DM(1.0), DM(2.0)), (DM(2.0), DM(2.0)), (DM(3.0), DM(2.0))],
        doc="Test less-than comparison"
    )

    # ======= FMIN test (element-wise min) =======
    x_fmin = MX.sym("x", 2, 2)
    y_fmin = MX.sym("y", 2, 2)
    Onnxtests.add_test(
        "fmin",
        Function("fmin_test", [x_fmin, y_fmin], [fmin(x_fmin, y_fmin)]),
        [(DM([[1, 5], [3, 2]]), DM([[4, 2], [1, 6]]))],
        doc="Test element-wise minimum"
    )

    # ======= FMAX test (element-wise max) =======
    x_fmax = MX.sym("x", 2, 2)
    y_fmax = MX.sym("y", 2, 2)
    Onnxtests.add_test(
        "fmax",
        Function("fmax_test", [x_fmax, y_fmax], [fmax(x_fmax, y_fmax)]),
        [(DM([[1, 5], [3, 2]]), DM([[4, 2], [1, 6]]))],
        doc="Test element-wise maximum"
    )

    # ======= Bilinear form test (x' * A * y) =======
    A_bl = MX.sym("A", 3, 3)
    x_bl = MX.sym("x", 3, 1)
    y_bl = MX.sym("y", 3, 1)
    Onnxtests.add_test(
        "bilin",
        Function("bilin_test", [A_bl, x_bl, y_bl], [bilin(A_bl, x_bl, y_bl)]),
        [(DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), DM([1, 2, 3]), DM([1, 0, 1]))],
        doc="Test bilinear form x'*A*y"
    )

    # ======= Rank-1 update test (A + alpha * x * y') =======
    A_r1 = MX.sym("A", 2, 2)
    alpha_r1 = MX.sym("alpha")
    x_r1 = MX.sym("x", 2, 1)
    y_r1 = MX.sym("y", 2, 1)
    Onnxtests.add_test(
        "rank1",
        Function("rank1_test", [A_r1, alpha_r1, x_r1, y_r1], [rank1(A_r1, alpha_r1, x_r1, y_r1)]),
        [(DM([[1, 2], [3, 4]]), DM(2.0), DM([1, 0]), DM([0, 1]))],
        doc="Test rank-1 update A + alpha*x*y'"
    )

    # ======= If-else-zero test (conditional) =======
    x_ie = MX.sym("x")
    y_ie = MX.sym("y")
    Onnxtests.add_test(
        "if_else_zero",
        Function("ifez_test", [x_ie, y_ie], [if_else(x_ie > 0, y_ie, 0)]),
        [(DM(1.0), DM(5.0)), (DM(-1.0), DM(5.0)), (DM(0.0), DM(3.0))],
        doc="Test if_else_zero conditional operation"
    )


# Register all tests
_register_op_tests()
_register_complex_tests()
_register_mx_ops_tests()
_register_tensor_op_tests()
_register_additional_tests()


if __name__ == '__main__':
    unittest.main()
