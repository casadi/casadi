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

# ============================================================================
# TEST DATA
# ============================================================================

# Unary operations: (name, expr_lambda, test_values, skip_onnxruntime)
# skip_onnxruntime=True for ops without float64 support in ONNX Runtime
UNARY_OPS = [
    ("sin", lambda x: sin(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)], False),
    ("cos", lambda x: cos(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)], True),   # No float64
    ("tan", lambda x: tan(x), [DM(0.0), DM(0.5), DM(1.0)], True),             # No float64
    ("exp", lambda x: exp(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)], False),
    ("log", lambda x: log(x), [DM(0.1), DM(1.0), DM(2.0), DM(10.0)], False),
    ("sqrt", lambda x: sqrt(x), [DM(0.0), DM(1.0), DM(4.0), DM(9.0)], False),
    ("asin", lambda x: asin(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)], True), # No float64
    ("acos", lambda x: acos(x), [DM(0.0), DM(0.5), DM(1.0)], True),           # No float64
    ("atan", lambda x: atan(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)], True), # No float64
    ("sinh", lambda x: sinh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)], True), # No float64
    ("cosh", lambda x: cosh(x), [DM(0.0), DM(0.5), DM(1.0)], True),           # No float64
    ("tanh", lambda x: tanh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)], False),
    ("asinh", lambda x: asinh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)], True), # No float64
    ("acosh", lambda x: acosh(x), [DM(1.0), DM(1.5), DM(2.0)], True),         # No float64
    ("atanh", lambda x: atanh(x), [DM(0.0), DM(0.5), DM(-0.5)], True),        # No float64
    ("ceil", lambda x: ceil(x), [DM(1.2), DM(1.5), DM(1.9), DM(-1.2)], False),
    ("floor", lambda x: floor(x), [DM(1.2), DM(1.5), DM(1.9), DM(-1.2)], False),
    ("fabs", lambda x: fabs(x), [DM(1.0), DM(-1.0), DM(0.0), DM(-5.5)], False),
    ("sign", lambda x: sign(x), [DM(1.0), DM(-1.0), DM(0.0), DM(5.5)], False),
    ("neg", lambda x: -x, [DM(1.0), DM(-1.0), DM(0.0), DM(5.5)], False),
    ("erf", lambda x: erf(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)], True),   # No float64
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
            t_export = Translator("onnx", {"verbose": False})
            t_export.load(casadi_func)
            t_export.save(onnx_file)

            # Import from ONNX
            t_import = Translator("onnx", {"verbose": False})
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
            t = Translator("onnx", {"verbose": False})
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
    def add_test(cls, test_name, func, test_inputs, doc=None, skip_onnxruntime=False):
        """
        Add a test (roundtrip + onnxruntime if available) to the test class.

        Args:
            test_name: Base name for the test (e.g., "unary_sin")
            func: CasADi Function to test
            test_inputs: List of test input values
            doc: Optional docstring for the test
            skip_onnxruntime: Skip onnxruntime validation (for ops without float64 support)
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

        # Create onnxruntime test (only if available and not skipped)
        if skip_onnxruntime:
            def skipped_test(self):
                self.skipTest("ONNX Runtime lacks float64 support for this op")
            skipped_test.__name__ = f"test_{test_name}_onnxruntime"
            skipped_test.__doc__ = doc or f"Test {test_name} via onnxruntime (skipped)"
            setattr(cls, skipped_test.__name__, skipped_test)
        elif HAS_ONNXRUNTIME:
            def make_onnxruntime_test(name, f, inputs):
                def test_method(self):
                    self.onnxruntime_test(name, f, inputs)
                return test_method

            onnxruntime_test = make_onnxruntime_test(test_name, func, test_inputs)
            onnxruntime_test.__name__ = f"test_{test_name}_onnxruntime"
            onnxruntime_test.__doc__ = doc or f"Test {test_name} via onnxruntime"
            setattr(cls, onnxruntime_test.__name__, onnxruntime_test)
        else:
            # If ONNX Runtime is not available, add a skipped test
            def skipped_test(self):
                self.skipTest("ONNX Runtime not available")

            skipped_test.__name__ = f"test_{test_name}_onnxruntime"
            skipped_test.__doc__ = doc or f"Test {test_name} via onnxruntime (skipped)"
            setattr(cls, skipped_test.__name__, skipped_test)


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
    for op_name, expr, test_values, skip_ort in UNARY_OPS:
        Onnxtests.add_test(
            f"unary_{op_name}",
            Function(f"test_{op_name}", [_x], [expr(_x)]),
            test_values,
            doc=f"Test unary {op_name} operation",
            skip_onnxruntime=skip_ort
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
        doc="Test function with no inputs",
        skip_onnxruntime=True  # Shape mismatch in ONNX Runtime
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

    # Control flow: if_else
    f_then = Function("then_branch", [_x], [_x + 1])
    f_else = Function("else_branch", [_x], [_x - 1])
    Onnxtests.add_test(
        "control_flow_if_else",
        Function.if_else("test_if", f_then, f_else).wrap(),
        [(DM(1), DM(5.0))],
        doc="Test if_else control flow",
        skip_onnxruntime=True  # SSA form issue in ONNX export
    )

    # Function hierarchy: simple
    # Note: skip_onnxruntime=True because FunctionProto export lacks opset imports
    f_inner = Function("inner", [_x], [_x + _x])
    Onnxtests.add_test(
        "function_hierarchy_simple",
        Function("outer", [_y], [f_inner(_y) + 1]).wrap(),
        [DM(5.0)],
        doc="Test nested function call: outer calls inner",
        skip_onnxruntime=True  # FunctionProto needs opset imports
    )

    # Function hierarchy: multiple calls
    f_double = Function("double", [_x], [_x + _x])
    Onnxtests.add_test(
        "function_hierarchy_multiple_calls",
        Function("outer_multi", [_y], [f_double(_y) + f_double(_y + 1)]).wrap(),
        [DM(3.0)],
        doc="Test outer function calling inner function twice",
        skip_onnxruntime=True  # FunctionProto needs opset imports
    )

    # Function hierarchy: deep (A calls B calls C)
    f_c = Function("func_c", [_x], [sin(_x)])
    f_b = Function("func_b", [_y], [f_c(_y) + 1])
    Onnxtests.add_test(
        "function_hierarchy_deep",
        Function("func_a", [_z], [f_b(_z) * 2]).wrap(),
        [DM(0.5)],
        doc="Test deep function hierarchy: A calls B calls C",
        skip_onnxruntime=True  # FunctionProto needs opset imports
    )

    # Control flow: map (Scan)
    f_map_body = Function("map_body", [_x], [sin(_x)])
    Onnxtests.add_test(
        "control_flow_map",
        f_map_body.map(3, "serial").wrap(),
        [DM([0.0, 1.0, 2.0])],
        doc="Test map (Scan) control flow",
        skip_onnxruntime=True  # Scan export issues
    )

    # Control flow: mapaccum (Loop)
    # TODO: mapaccum outputs all intermediate values to a single array, but ONNX export
    # currently only keeps the last value. Needs Concat support to collect all outputs.
    # f_accum_body = Function("accum_body", [_x], [_x + 1])
    # Onnxtests.add_test(
    #     "control_flow_mapaccum",
    #     f_accum_body.mapaccum("test_accum", 3).wrap(),
    #     [DM(0.0)],
    #     doc="Test mapaccum (Loop) control flow"
    # )


# Register all tests
_register_op_tests()
_register_complex_tests()


if __name__ == '__main__':
    unittest.main()
