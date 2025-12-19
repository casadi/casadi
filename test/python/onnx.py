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
# TEST CASE REGISTRY (Data-Driven Approach)
# ============================================================================
# To add a new operation test, just add one line to the appropriate list!

UNARY_OPS = [
    # (name, expression_builder, test_input_values)
    ("sin", lambda x: sin(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("cos", lambda x: cos(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("tan", lambda x: tan(x), [DM(0.0), DM(0.5), DM(1.0)]),
    ("exp", lambda x: exp(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)]),
    ("log", lambda x: log(x), [DM(0.1), DM(1.0), DM(2.0), DM(10.0)]),  # x > 0
    ("sqrt", lambda x: sqrt(x), [DM(0.0), DM(1.0), DM(4.0), DM(9.0)]),
    ("asin", lambda x: asin(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("acos", lambda x: acos(x), [DM(0.0), DM(0.5), DM(1.0)]),
    ("atan", lambda x: atan(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)]),
    ("sinh", lambda x: sinh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("cosh", lambda x: cosh(x), [DM(0.0), DM(0.5), DM(1.0)]),
    ("tanh", lambda x: tanh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-1.0)]),
    ("asinh", lambda x: asinh(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
    ("acosh", lambda x: acosh(x), [DM(1.0), DM(1.5), DM(2.0)]),  # x >= 1
    ("atanh", lambda x: atanh(x), [DM(0.0), DM(0.5), DM(-0.5)]),  # |x| < 1
    ("ceil", lambda x: ceil(x), [DM(1.2), DM(1.5), DM(1.9), DM(-1.2)]),
    ("floor", lambda x: floor(x), [DM(1.2), DM(1.5), DM(1.9), DM(-1.2)]),
    ("fabs", lambda x: fabs(x), [DM(1.0), DM(-1.0), DM(0.0), DM(-5.5)]),
    ("sign", lambda x: sign(x), [DM(1.0), DM(-1.0), DM(0.0), DM(5.5)]),
    ("neg", lambda x: -x, [DM(1.0), DM(-1.0), DM(0.0), DM(5.5)]),
    ("erf", lambda x: erf(x), [DM(0.0), DM(0.5), DM(1.0), DM(-0.5)]),
]

BINARY_OPS = [
    # (name, expression_builder, test_input_pairs)
    ("add", lambda x, y: x + y, [(DM(2.0), DM(3.0)), (DM(-1.0), DM(4.0)), (DM(0.0), DM(0.0))]),
    ("sub", lambda x, y: x - y, [(DM(5.5), DM(2.2)), (DM(10.0), DM(3.0)), (DM(0.0), DM(0.0))]),
    ("mul", lambda x, y: x * y, [(DM(2.0), DM(3.0)), (DM(-2.0), DM(4.0)), (DM(0.0), DM(5.0))]),
    ("div", lambda x, y: x / y, [(DM(10.0), DM(2.0)), (DM(7.0), DM(3.0)), (DM(1.0), DM(1.0))]),
    ("pow", lambda x, y: x ** y, [(DM(2.0), DM(3.0)), (DM(4.0), DM(0.5)), (DM(10.0), DM(2.0))]),
]

MATRIX_OPS = [
    # (name, expression_builder, input_shapes, test_value_generator)
    (
        "matmul_2x3_3x2",
        lambda A, B: mtimes(A, B),
        [(2, 3), (3, 2)],
        lambda: (DM([[1, 2, 3], [4, 5, 6]]), DM([[1, 2], [3, 4], [5, 6]])),
    ),
    (
        "matmul_3x1_1x3",
        lambda A, B: mtimes(A, B),
        [(3, 1), (1, 3)],
        lambda: (DM([[1], [2], [3]]), DM([[4, 5, 6]])),
    ),
    (
        "transpose_2x3",
        lambda A: A.T,
        [(2, 3)],
        lambda: (DM([[1, 2, 3], [4, 5, 6]]),),
    ),
]

SPECIAL_OPS = [
    # Operations requiring unique test logic (not in parametrized tests)
    "identity",
    "constant",
    "concat",
    "split",
]

ERROR_CASES = [
    # (name, expression_builder, expected_error_substring)
    ("unsupported_op", lambda x: x, "unsupported operation"),  # Placeholder for actual unsupported op
]


# ============================================================================
# HELPER METHODS (Reusable Test Logic)
# ============================================================================

class Onnxtests(casadiTestCase):
    """ONNX translator tests"""

    def roundtrip_test(self, op_name, casadi_func, test_inputs):
        """
        Execute full roundtrip test: CasADi → ONNX → CasADi → validate

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

                self.checkarray(result_original, result_imported, digits=10)

            # Also validate with ONNX Runtime if available
            first_input = test_inputs[0] if not isinstance(test_inputs[0], tuple) else test_inputs[0]
            self.validate_with_onnxruntime(onnx_file, casadi_func, first_input)

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

            # Validate with ONNX Runtime
            test_val = test_inputs[0] if not isinstance(test_inputs[0], tuple) else test_inputs[0]
            self.validate_with_onnxruntime(onnx_file, casadi_func, test_val)

        finally:
            if os.path.exists(onnx_file):
                os.remove(onnx_file)

    def validate_with_onnxruntime(self, onnx_file, casadi_function, test_inputs):
        """
        Validate CasADi function output against ONNX Runtime
        (Existing method from onnx.py - kept as-is for compatibility)
        """
        try:
            import onnxruntime as ort
        except ImportError:
            # Graceful degradation - skip if onnxruntime not available
            return

        # Prepare inputs
        if not isinstance(test_inputs, (list, tuple)):
            test_inputs = [test_inputs]

        # Convert to numpy format for ONNX Runtime
        import numpy as np
        onnx_inputs = {}
        for i, inp in enumerate(test_inputs):
            input_name = casadi_function.name_in(i)
            if input_name == "":
                input_name = f"input_{i}"
            onnx_inputs[input_name] = np.array(inp).astype(np.float64)

        # Run CasADi
        casadi_outputs = casadi_function(*test_inputs)
        if not isinstance(casadi_outputs, (list, tuple)):
            casadi_outputs = [casadi_outputs]

        # Run ONNX Runtime
        session = ort.InferenceSession(onnx_file)
        onnx_outputs = session.run(None, onnx_inputs)

        # Compare
        for i, (onnx_out, casadi_out) in enumerate(zip(onnx_outputs, casadi_outputs)):
            self.checkarray(DM(onnx_out), casadi_out, digits=10)


    # ========================================================================
    # PARAMETRIZED TESTS (Auto-generated from registry)
    # ========================================================================
    # Note: Individual test methods are added dynamically below using
    # _add_parametrized_tests(). Each operation gets its own test method.
    # E.g., test_unary_roundtrip_sin(), test_unary_roundtrip_cos(), etc.
    #
    # This is pure unittest - no external libraries needed!


    # ========================================================================
    # SPECIAL TESTS (Unique Logic, Not Data-Driven)
    # ========================================================================

    def test_complex_expression_graph(self):
        """
        Test complex expression combining multiple operations

        Example of a special test that doesn't fit the data-driven pattern
        """
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        # Create complex expression
        x = MX.sym("x")
        y = MX.sym("y")

        # Complex expression: (sin(x) + cos(y)) * exp(x - y)
        expr = (sin(x) + cos(y)) * exp(x - y)
        f = Function("complex_expr", [x, y], [expr])

        # Test values
        test_values = [
            (DM(0.5), DM(1.0)),
            (DM(1.0), DM(0.5)),
            (DM(0.0), DM(0.0)),
        ]

        # Run roundtrip test
        self.roundtrip_test("complex_expr", f, test_values)

    def test_multiple_outputs(self):
        """Test function with multiple outputs"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        x = MX.sym("x")
        y = MX.sym("y")

        # Function with 3 outputs
        f = Function("multi_out", [x, y], [x + y, x - y, x * y])

        test_values = [
            (DM(3.0), DM(2.0)),
            (DM(5.0), DM(1.0)),
        ]

        self.roundtrip_test("multi_out", f, test_values)

    def test_different_tensor_shapes(self):
        """Test operations with different tensor shapes"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        shapes_to_test = [
            # (shape, test_value)
            ((1, 1), DM([[2.0]])),      # Scalar (1x1)
            ((3, 1), DM([[1.0], [2.0], [3.0]])),  # Column vector
            ((1, 3), DM([[1.0, 2.0, 3.0]])),      # Row vector
            ((2, 2), DM([[1.0, 2.0], [3.0, 4.0]])),  # 2x2 matrix
        ]

        for shape, test_val in shapes_to_test:
            with self.subTest(shape=f"{shape[0]}x{shape[1]}"):
                x = MX.sym("x", shape[0], shape[1])
                f = Function(f"sin_{shape[0]}x{shape[1]}", [x], [sin(x)])

                self.roundtrip_test(f"sin_{shape[0]}x{shape[1]}", f, [test_val])

    def test_empty_function(self):
        """Test function with no inputs/outputs"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        # Create empty function (just returns a constant)
        # Use MX to ensure MXFunction (not SXFunction)
        f = Function("empty", [], [MX(42.0)])

        onnx_file = "test_empty.onnx"
        try:
            t = translator("onnx", {"verbose": False})
            t.load(f)
            t.save(onnx_file)

            t2 = translator("onnx", {"verbose": False})
            t2.load(onnx_file)
            f2 = t2.create("empty_imported")

            # Validate
            self.checkarray(f(), f2(), digits=10)
        finally:
            if os.path.exists(onnx_file):
                os.remove(onnx_file)

    # ========================================================================
    # Control Flow Tests
    # ========================================================================

    def test_control_flow_if_else(self):
        """Test if_else control flow roundtrip"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        # Create a simple if_else function: if x > 0: x+1, else: x-1
        x = MX.sym("x")
        f_then = Function("then_branch", [x], [x + 1])
        f_else = Function("else_branch", [x], [x - 1])
        f_if = Function.if_else("test_if", f_then, f_else)

        # Wrap to MXFunction for ONNX export
        f_if = f_if.wrap()

        # Test with positive value (should use then branch: 5 + 1 = 6)
        self.roundtrip_test("if_else", f_if, [DM(5.0)])

    def test_control_flow_map(self):
        """Test map (Scan) control flow roundtrip"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        # Create a simple map function: apply sin to each element
        x = MX.sym("x")
        f_base = Function("map_body", [x], [sin(x)])
        f_map = f_base.map(3, "serial")

        # Wrap to MXFunction for ONNX export
        f_map = f_map.wrap()

        # Test with array of 3 elements
        self.roundtrip_test("map", f_map, [DM([0.0, 1.0, 2.0])])

    def test_control_flow_mapaccum(self):
        """Test mapaccum (Loop) control flow roundtrip"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        # Create a simple mapaccum function: accumulate by adding 1
        x = MX.sym("x")
        f_base = Function("accum_body", [x], [x + 1])
        f_accum = f_base.mapaccum("test_accum", 3)

        # Wrap to MXFunction for ONNX export
        f_accum = f_accum.wrap()

        # Test with initial value 0 (should produce: 1, 2, 3)
        self.roundtrip_test("mapaccum", f_accum, [DM(0.0)])

    def test_conditional_loop_termination(self):
        """Test Loop with conditional termination (import from ONNX)"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        try:
            import onnx
            from onnx import helper, TensorProto, numpy_helper
        except ImportError:
            self.skipTest("ONNX Python package not available")

        # Create an ONNX Loop that counts from 0 until reaching a threshold
        # Loop body: (iter, cond, counter) -> (cond_out, counter+1)
        #   cond_out = (counter < threshold)

        # Create the loop body subgraph
        # Inputs: iter_num, cond_in, counter
        # Outputs: cond_out, counter_updated

        body_iter = helper.make_tensor_value_info("iter", TensorProto.DOUBLE, [1, 1])
        body_cond_in = helper.make_tensor_value_info("cond_in", TensorProto.DOUBLE, [1, 1])
        body_counter = helper.make_tensor_value_info("counter", TensorProto.DOUBLE, [1, 1])

        # Constant: threshold = 5
        threshold = helper.make_tensor("threshold", TensorProto.DOUBLE, [1, 1], [5.0])

        # Constant: one = 1
        one = helper.make_tensor("one", TensorProto.DOUBLE, [1, 1], [1.0])

        # counter_updated = counter + 1
        add_node = helper.make_node("Add", ["counter", "one"], ["counter_updated"])

        # cond_out = counter < threshold
        # Note: ONNX Less returns bool, but we need double for consistency
        # We'll check counter_updated < threshold (so loop runs while counter < 5)
        less_node = helper.make_node("Less", ["counter_updated", "threshold"], ["cond_out_bool"])

        # Convert bool to double (Cast node)
        cast_node = helper.make_node("Cast", ["cond_out_bool"], ["cond_out"],
                                     to=TensorProto.DOUBLE)

        body_cond_out = helper.make_tensor_value_info("cond_out", TensorProto.DOUBLE, [1, 1])
        body_counter_out = helper.make_tensor_value_info("counter_updated", TensorProto.DOUBLE, [1, 1])

        # Create the body graph
        body_graph = helper.make_graph(
            [add_node, less_node, cast_node],
            "loop_body",
            [body_iter, body_cond_in, body_counter],
            [body_cond_out, body_counter_out],
            [threshold, one]
        )

        # Create the main graph with Loop operator
        # Inputs: initial counter value
        graph_input = helper.make_tensor_value_info("initial_counter", TensorProto.DOUBLE, [1, 1])

        # Loop inputs: max_iter (10), initial_cond (true), counter (0)
        max_iter = helper.make_tensor("max_iter", TensorProto.INT64, [], [10])
        init_cond = helper.make_tensor("init_cond", TensorProto.DOUBLE, [1, 1], [1.0])  # true

        # Loop node
        loop_node = helper.make_node(
            "Loop",
            ["max_iter", "init_cond", "initial_counter"],
            ["final_counter"],
            body=body_graph
        )

        # Output
        graph_output = helper.make_tensor_value_info("final_counter", TensorProto.DOUBLE, [1, 1])

        # Create the graph
        graph = helper.make_graph(
            [loop_node],
            "conditional_loop_test",
            [graph_input],
            [graph_output],
            [max_iter, init_cond]
        )

        # Create the model
        model = helper.make_model(graph, producer_name="casadi_test")

        # Save to temp file
        import tempfile
        import os
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.onnx', delete=False) as f:
            onnx.save(model, f.name)
            temp_file = f.name

        try:
            # Import the ONNX model
            t = Translator.load("onnx", temp_file)
            f = t.create("conditional_loop")

            # Test: starting from 0, should count to 5 (when counter becomes 5, condition becomes false)
            result = f(DM(0.0))
            expected = DM(5.0)

            # Verify the result
            self.assertTrue(result.is_dense(), "Result should be dense")
            self.assertEqual(result.size1(), 1, "Result should be 1x1")
            self.assertEqual(result.size2(), 1, "Result should be 1x1")
            self.assertAlmostEqual(float(result), float(expected), places=5,
                                 msg=f"Loop should terminate at threshold=5, got {float(result)}")

        finally:
            # Clean up temp file
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_import_float32_constant(self):
        """Test importing FLOAT (32-bit) tensors from ONNX"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        try:
            import onnx
            from onnx import helper, TensorProto
        except ImportError:
            self.skipTest("ONNX Python package not available")

        # Create ONNX model with FLOAT constant
        import tempfile
        import os

        float_values = [1.5, 2.5, 3.5, 4.5]
        constant_node = helper.make_node(
            'Constant',
            inputs=[],
            outputs=['output'],
            value=helper.make_tensor(
                name='const_tensor',
                data_type=TensorProto.FLOAT,
                dims=[2, 2],
                vals=float_values
            )
        )

        graph_output = helper.make_tensor_value_info('output', TensorProto.FLOAT, [2, 2])

        graph = helper.make_graph([constant_node], 'float_test', [], [graph_output])
        model = helper.make_model(graph, producer_name='casadi_test')

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.onnx', delete=False) as f:
            onnx.save(model, f.name)
            temp_file = f.name

        try:
            t = Translator.load("onnx", temp_file)
            func = t.create("float_import_test")
            result = func()

            # Verify values (should be promoted to double)
            self.assertEqual(result.size1(), 2)
            self.assertEqual(result.size2(), 2)
            for i in range(4):
                self.assertAlmostEqual(float(result[i]), float_values[i], places=6)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_import_int32_constant(self):
        """Test importing INT32 tensors from ONNX"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        try:
            import onnx
            from onnx import helper, TensorProto
        except ImportError:
            self.skipTest("ONNX Python package not available")

        import tempfile
        import os

        int_values = [10, 20, 30, 40, 50, 60]
        constant_node = helper.make_node(
            'Constant',
            inputs=[],
            outputs=['output'],
            value=helper.make_tensor(
                name='const_tensor',
                data_type=TensorProto.INT32,
                dims=[2, 3],
                vals=int_values
            )
        )

        graph_output = helper.make_tensor_value_info('output', TensorProto.INT32, [2, 3])

        graph = helper.make_graph([constant_node], 'int32_test', [], [graph_output])
        model = helper.make_model(graph, producer_name='casadi_test')

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.onnx', delete=False) as f:
            onnx.save(model, f.name)
            temp_file = f.name

        try:
            t = Translator.load("onnx", temp_file)
            func = t.create("int32_import_test")
            result = func()

            # Verify values (should be converted to double)
            self.assertEqual(result.size1(), 2)
            self.assertEqual(result.size2(), 3)
            for i in range(6):
                self.assertAlmostEqual(float(result[i]), float(int_values[i]), places=1)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_import_int64_constant(self):
        """Test importing INT64 tensors from ONNX"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        try:
            import onnx
            from onnx import helper, TensorProto
        except ImportError:
            self.skipTest("ONNX Python package not available")

        import tempfile
        import os

        # Use small INT64 values (within double precision range)
        int64_values = [100, 200, 300, 400]
        constant_node = helper.make_node(
            'Constant',
            inputs=[],
            outputs=['output'],
            value=helper.make_tensor(
                name='const_tensor',
                data_type=TensorProto.INT64,
                dims=[1, 4],
                vals=int64_values
            )
        )

        graph_output = helper.make_tensor_value_info('output', TensorProto.INT64, [1, 4])

        graph = helper.make_graph([constant_node], 'int64_test', [], [graph_output])
        model = helper.make_model(graph, producer_name='casadi_test')

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.onnx', delete=False) as f:
            onnx.save(model, f.name)
            temp_file = f.name

        try:
            t = Translator.load("onnx", temp_file)
            func = t.create("int64_import_test")
            result = func()

            # Verify values (should be converted to double)
            self.assertEqual(result.size1(), 1)
            self.assertEqual(result.size2(), 4)
            for i in range(4):
                self.assertAlmostEqual(float(result[i]), float(int64_values[i]), places=1)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_import_bool_constant(self):
        """Test importing BOOL tensors from ONNX"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        try:
            import onnx
            from onnx import helper, TensorProto
        except ImportError:
            self.skipTest("ONNX Python package not available")

        import tempfile
        import os

        bool_values = [True, False, True, False, True, True]
        constant_node = helper.make_node(
            'Constant',
            inputs=[],
            outputs=['output'],
            value=helper.make_tensor(
                name='const_tensor',
                data_type=TensorProto.BOOL,
                dims=[3, 2],
                vals=bool_values
            )
        )

        graph_output = helper.make_tensor_value_info('output', TensorProto.BOOL, [3, 2])

        graph = helper.make_graph([constant_node], 'bool_test', [], [graph_output])
        model = helper.make_model(graph, producer_name='casadi_test')

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.onnx', delete=False) as f:
            onnx.save(model, f.name)
            temp_file = f.name

        try:
            t = Translator.load("onnx", temp_file)
            func = t.create("bool_import_test")
            result = func()

            # Verify values (should be converted to 0.0/1.0)
            self.assertEqual(result.size1(), 3)
            self.assertEqual(result.size2(), 2)
            expected = [1.0, 0.0, 1.0, 0.0, 1.0, 1.0]
            for i in range(6):
                self.assertAlmostEqual(float(result[i]), expected[i], places=1)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_import_mixed_types(self):
        """Test ONNX model with multiple tensor types"""
        if not has_translator("onnx"):
            self.skipTest("ONNX translator not available")

        try:
            import onnx
            from onnx import helper, TensorProto
        except ImportError:
            self.skipTest("ONNX Python package not available")

        import tempfile
        import os

        # Create nodes with different types
        # FLOAT constant
        float_const = helper.make_node(
            'Constant',
            inputs=[],
            outputs=['float_val'],
            value=helper.make_tensor('float_tensor', TensorProto.FLOAT, [1, 1], [2.5])
        )

        # INT32 constant
        int_const = helper.make_node(
            'Constant',
            inputs=[],
            outputs=['int_val'],
            value=helper.make_tensor('int_tensor', TensorProto.INT32, [1, 1], [3])
        )

        # Add them (both will be double after import)
        add_node = helper.make_node('Add', ['float_val', 'int_val'], ['output'])

        graph_output = helper.make_tensor_value_info('output', TensorProto.DOUBLE, [1, 1])

        graph = helper.make_graph(
            [float_const, int_const, add_node],
            'mixed_types_test',
            [],
            [graph_output]
        )
        model = helper.make_model(graph, producer_name='casadi_test')

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.onnx', delete=False) as f:
            onnx.save(model, f.name)
            temp_file = f.name

        try:
            t = Translator.load("onnx", temp_file)
            func = t.create("mixed_types_test")
            result = func()

            # Verify result: 2.5 + 3 = 5.5
            self.assertAlmostEqual(float(result), 5.5, places=5)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


# ============================================================================
# DYNAMIC TEST GENERATION (Parametrized Tests for unittest)
# ============================================================================

def _add_parametrized_tests():
    """
    Dynamically generate individual test methods for each operation.

    This function creates test methods like:
    - test_unary_roundtrip_sin()
    - test_unary_roundtrip_cos()
    - test_binary_roundtrip_add()
    - etc.

    Pure unittest - no external libraries needed!
    """

    def make_unary_test(op_name, expr, test_values):
        """Factory function to create a unary operation test"""
        def test(self):
            if not has_translator("onnx"):
                self.skipTest("ONNX translator not available")

            x = MX.sym("x")
            f = Function(f"test_{op_name}", [x], [expr(x)])
            self.roundtrip_test(op_name, f, test_values)

        test.__name__ = f"test_unary_roundtrip_{op_name}"
        test.__doc__ = f"Test {op_name} operation via roundtrip"
        return test

    def make_binary_test(op_name, expr, test_values):
        """Factory function to create a binary operation test"""
        def test(self):
            if not has_translator("onnx"):
                self.skipTest("ONNX translator not available")

            x = MX.sym("x")
            y = MX.sym("y")
            f = Function(f"test_{op_name}", [x, y], [expr(x, y)])
            self.roundtrip_test(op_name, f, test_values)

        test.__name__ = f"test_binary_roundtrip_{op_name}"
        test.__doc__ = f"Test {op_name} operation via roundtrip"
        return test

    def make_matrix_test(op_name, expr, shapes, value_gen):
        """Factory function to create a matrix operation test"""
        def test(self):
            if not has_translator("onnx"):
                self.skipTest("ONNX translator not available")

            if len(shapes) == 1:
                A = MX.sym("A", shapes[0][0], shapes[0][1])
                f = Function(f"test_{op_name}", [A], [expr(A)])
            else:
                A = MX.sym("A", shapes[0][0], shapes[0][1])
                B = MX.sym("B", shapes[1][0], shapes[1][1])
                f = Function(f"test_{op_name}", [A, B], [expr(A, B)])

            test_values = [value_gen()]
            self.roundtrip_test(op_name, f, test_values)

        test.__name__ = f"test_matrix_roundtrip_{op_name}"
        test.__doc__ = f"Test {op_name} matrix operation via roundtrip"
        return test

    # Add unary operation tests
    for op_name, expr, test_values in UNARY_OPS:
        test_method = make_unary_test(op_name, expr, test_values)
        setattr(Onnxtests, test_method.__name__, test_method)

    # Add binary operation tests
    for op_name, expr, test_values in BINARY_OPS:
        test_method = make_binary_test(op_name, expr, test_values)
        setattr(Onnxtests, test_method.__name__, test_method)

    # Add matrix operation tests
    for op_name, expr, shapes, value_gen in MATRIX_OPS:
        test_method = make_matrix_test(op_name, expr, shapes, value_gen)
        setattr(Onnxtests, test_method.__name__, test_method)

# Generate all parametrized tests
_add_parametrized_tests()


if __name__ == '__main__':
    unittest.main()
