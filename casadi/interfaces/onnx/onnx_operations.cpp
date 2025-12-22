/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "onnx_translator.hpp"

/// \cond INTERNAL
namespace casadi {

  MX OnnxTranslator::process_node_operation(
      const std::string& op_type,
      const onnx::NodeProto& node,
      const std::vector<MX>& node_inputs) {

    MX output;

    // Try simple operations using centralized lookup table
    const OpMapping* mapping = get_op_mapping_by_name(op_type);
    if (mapping) {
      const MX& x = node_inputs[0];
      if (mapping->arity == 1) {
        casadi_assert(node_inputs.size() >= 1, op_type + " requires 1 input");
        switch (mapping->casadi_op) {
          case OP_SIN: return sin(x);
          case OP_COS: return cos(x);
          case OP_TAN: return tan(x);
          case OP_ASIN: return asin(x);
          case OP_ACOS: return acos(x);
          case OP_ATAN: return atan(x);
          case OP_SINH: return sinh(x);
          case OP_COSH: return cosh(x);
          case OP_TANH: return tanh(x);
          case OP_ASINH: return asinh(x);
          case OP_ACOSH: return acosh(x);
          case OP_ATANH: return atanh(x);
          case OP_EXP: return exp(x);
          case OP_LOG: return log(x);
          case OP_SQRT: return sqrt(x);
          case OP_NEG: return -x;
          case OP_FABS: return fabs(x);
          case OP_CEIL: return ceil(x);
          case OP_FLOOR: return floor(x);
          case OP_SIGN: return sign(x);
          case OP_ERF: return erf(x);
          case OP_INV: return 1.0 / x;
          case OP_NOT: return !x;
          case OP_TRANSPOSE: return x.T();
          case OP_NORM1: return norm_1(x);
          case OP_NORM2: return norm_2(x);
          case OP_NORMF: return norm_fro(x);
          case OP_MMIN: return mmin(x);
          case OP_MMAX: return mmax(x);
          default: break;
        }
      } else if (mapping->arity == 2) {
        casadi_assert(node_inputs.size() >= 2, op_type + " requires 2 inputs");
        const MX& y = node_inputs[1];
        switch (mapping->casadi_op) {
          case OP_ADD: return x + y;
          case OP_SUB: return x - y;
          case OP_MUL: return x * y;
          case OP_DIV: return x / y;
          case OP_POW: return pow(x, y);
          case OP_OR: return logic_or(x, y);
          case OP_AND: return logic_and(x, y);
          default: break;
        }
      }
    }

    // Operations not in the simple mapping tables
    // These require special handling (attributes, multiple inputs, etc.)

    if (op_type == "Less") {
      casadi_assert(node_inputs.size() >= 2, "Less requires 2 inputs");
      // Less: returns 1.0 if input[0] < input[1], 0.0 otherwise
      // Use CasADi's if_else: if (a < b) then 1.0 else 0.0
      output = if_else(node_inputs[0] < node_inputs[1], MX(1.0), MX(0.0));

    // Utilities
    } else if (op_type == "Identity") {
      casadi_assert(node_inputs.size() >= 1, "Identity requires 1 input");
      output = node_inputs[0];

    } else if (op_type == "Cast") {
      // Cast operation: type conversion
      // We support importing DOUBLE, FLOAT, INT32, INT64, BOOL from ONNX,
      // but all types are converted to double at the tensor_to_dm() boundary.
      // Within CasADi's symbolic framework, everything is treated as double.
      // Therefore, Cast is effectively an identity operation during symbolic computation.
      casadi_assert(node_inputs.size() >= 1, "Cast requires 1 input");
      output = node_inputs[0];
      // Note: ONNX Cast has a 'to' attribute specifying target type.
      // Type conversions happen at ONNX import/export boundaries, not during symbolic ops.

    // Tensor operations
    } else if (op_type == "MatMul") {
      casadi_assert(node_inputs.size() >= 2, "MatMul requires 2 inputs");
      output = mtimes(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Constant") {
      // Extract value attribute
      const onnx::AttributeProto* value_attr = nullptr;
      for (int a = 0; a < node.attribute_size(); ++a) {
        if (node.attribute(a).name() == "value") {
          value_attr = &node.attribute(a);
          break;
        }
      }
      casadi_assert(value_attr != nullptr,
                   "Constant node must have 'value' attribute");
      const onnx::TensorProto& tensor = value_attr->t();
      DM dm_const = tensor_to_dm(tensor);
      output = MX(dm_const);

    // Complex tensor operations
    } else if (op_type == "Transpose") {
      casadi_assert(node_inputs.size() >= 1,
                    "Transpose operation requires 1 input");
      output = node_inputs[0].T();

    } else if (op_type == "Reshape") {
      casadi_assert(node_inputs.size() >= 2,
                    "Reshape operation requires 2 inputs (data and shape)");
      // Second input is the target shape - should be a constant
      casadi_assert(node_inputs[1].is_constant(),
                    "Reshape shape must be a constant");
      DM shape_dm = static_cast<DM>(node_inputs[1]);
      // Extract dimensions (assuming 2D for now)
      casadi_int new_rows = static_cast<casadi_int>(shape_dm(0).scalar());
      casadi_int new_cols = (shape_dm.numel() > 1) ?
                             static_cast<casadi_int>(shape_dm(1).scalar()) : 1;
      output = reshape(node_inputs[0], new_rows, new_cols);

    } else if (op_type == "Concat") {
      // Get axis attribute
      casadi_int axis = 0;
      for (int a = 0; a < node.attribute_size(); ++a) {
        if (node.attribute(a).name() == "axis") {
          axis = node.attribute(a).i();
          break;
        }
      }

      // Convert node_inputs vector to inputs for concat
      if (axis == 0) {
        // Vertical concatenation
        output = vertcat(node_inputs);
      } else if (axis == 1) {
        // Horizontal concatenation
        output = horzcat(node_inputs);
      } else {
        casadi_error("Concat with axis=" + std::to_string(axis) +
                     " not supported. Only axis=0 (vertcat) and axis=1 (horzcat) are supported.");
      }

    } else if (op_type == "Slice") {
      // Slice operation - extract sub-tensor
      // ONNX Slice has inputs: data, starts, ends, [axes], [steps]
      casadi_assert(node_inputs.size() >= 3, "Slice requires at least 3 inputs (data, starts, ends)");

      MX data = node_inputs[0];

      // Extract starts and ends (should be constants)
      casadi_assert(node_inputs[1].is_constant() && node_inputs[2].is_constant(),
                    "Slice starts and ends must be constants");

      DM starts_dm = static_cast<DM>(node_inputs[1]);
      DM ends_dm = static_cast<DM>(node_inputs[2]);

      // Extract axes if provided (default: [0, 1, ...])
      std::vector<casadi_int> axes;
      if (node_inputs.size() >= 4 && !node_inputs[3].is_empty()) {
        casadi_assert(node_inputs[3].is_constant(), "Slice axes must be constant");
        DM axes_dm = static_cast<DM>(node_inputs[3]);
        for (casadi_int k = 0; k < axes_dm.numel(); ++k) {
          axes.push_back(static_cast<casadi_int>(axes_dm(k).scalar()));
        }
      } else {
        // Default axes
        for (casadi_int k = 0; k < starts_dm.numel(); ++k) {
          axes.push_back(k);
        }
      }

      // Extract steps if provided (default: all 1s)
      std::vector<casadi_int> steps;
      if (node_inputs.size() >= 5 && !node_inputs[4].is_empty()) {
        casadi_assert(node_inputs[4].is_constant(), "Slice steps must be constant");
        DM steps_dm = static_cast<DM>(node_inputs[4]);
        for (casadi_int k = 0; k < steps_dm.numel(); ++k) {
          steps.push_back(static_cast<casadi_int>(steps_dm(k).scalar()));
        }
      } else {
        // Default steps: all 1s
        for (casadi_int k = 0; k < starts_dm.numel(); ++k) {
          steps.push_back(1);
        }
      }

      // Simplified implementation: only handle basic 2D slicing with step=1
      // Full implementation would need to handle arbitrary dimensions and steps
      casadi_assert(axes.size() <= 2, "Slice: only 2D slicing supported for now");
      casadi_assert(steps[0] == 1, "Slice: only step=1 supported for now");

      casadi_int start0 = static_cast<casadi_int>(starts_dm(0).scalar());
      casadi_int end0 = static_cast<casadi_int>(ends_dm(0).scalar());

      if (axes.size() == 1) {
        // Single axis slice
        if (axes[0] == 0) {
          // Row slice
          Slice row_slice(start0, end0);
          output = data(row_slice, Slice());
        } else {
          // Column slice
          Slice col_slice(start0, end0);
          output = data(Slice(), col_slice);
        }
      } else {
        // Two axis slice
        casadi_int start1 = static_cast<casadi_int>(starts_dm(1).scalar());
        casadi_int end1 = static_cast<casadi_int>(ends_dm(1).scalar());

        Slice row_slice(start0, end0);
        Slice col_slice(start1, end1);
        output = data(row_slice, col_slice);
      }

    } else if (op_type == "GatherElements") {
      // GatherElements - advanced indexing operation
      // Requires data tensor and indices tensor
      casadi_assert(node_inputs.size() >= 2, "GatherElements requires data and indices");

      // Get axis attribute
      int axis = 0;
      for (int a = 0; a < node.attribute_size(); ++a) {
        if (node.attribute(a).name() == "axis") {
          axis = node.attribute(a).i();
          break;
        }
      }

      // Full implementation requires advanced CasADi indexing capabilities
      casadi_error("ONNX import: GatherElements operation is not yet supported. "
                  "Cannot import ONNX models containing GatherElements nodes.");

    } else {
      casadi_error("Unsupported operation '" + op_type + "'");
    }

    return output;
  }

  // ========== Operation Mapping ==========
  //
  // Centralized mapping between CasADi opcodes and ONNX operation names.
  // Used by both export (CasADi → ONNX) and import (ONNX → CasADi).
  //
  // Single source of truth for all simple operations that have direct
  // CasADi ↔ ONNX equivalents.
  //
  // OpMapping struct is defined in onnx_translator.hpp

  static const OpMapping op_map[] = {
    // Unary operations
    {OP_SIN, "Sin", 1},
    {OP_COS, "Cos", 1},
    {OP_TAN, "Tan", 1},
    {OP_ASIN, "Asin", 1},
    {OP_ACOS, "Acos", 1},
    {OP_ATAN, "Atan", 1},
    {OP_SINH, "Sinh", 1},
    {OP_COSH, "Cosh", 1},
    {OP_TANH, "Tanh", 1},
    {OP_ASINH, "Asinh", 1},
    {OP_ACOSH, "Acosh", 1},
    {OP_ATANH, "Atanh", 1},
    {OP_EXP, "Exp", 1},
    {OP_LOG, "Log", 1},
    {OP_SQRT, "Sqrt", 1},
    {OP_NEG, "Neg", 1},
    {OP_FABS, "Abs", 1},
    {OP_CEIL, "Ceil", 1},
    {OP_FLOOR, "Floor", 1},
    {OP_SIGN, "Sign", 1},
    {OP_ERF, "Erf", 1},
    {OP_INV, "Reciprocal", 1},
    {OP_NOT, "Not", 1},
    {OP_TRANSPOSE, "Transpose", 1},
    {OP_NORM1, "ReduceL1", 1},
    {OP_NORM2, "ReduceL2", 1},
    {OP_NORMF, "ReduceL2", 1},
    {OP_MMIN, "ReduceMin", 1},
    {OP_MMAX, "ReduceMax", 1},
    // Binary operations
    {OP_ADD, "Add", 2},
    {OP_SUB, "Sub", 2},
    {OP_MUL, "Mul", 2},
    {OP_DIV, "Div", 2},
    {OP_POW, "Pow", 2},
    {OP_OR, "Or", 2},
    {OP_AND, "And", 2},
  };

  // Lookup by CasADi opcode (for export)
  const OpMapping* get_op_mapping(casadi_int op) {
    for (const auto& entry : op_map) {
      if (entry.casadi_op == op) return &entry;
    }
    return nullptr;
  }

  // Lookup by ONNX name (for import)
  const OpMapping* get_op_mapping_by_name(const std::string& onnx_name) {
    for (const auto& entry : op_map) {
      if (entry.onnx_name == onnx_name) return &entry;
    }
    return nullptr;
  }

  // ========== Export Helpers ==========

  // Callback-based implementations (main logic)
  onnx::NodeProto* create_binary_node(
      AddNodeFn add_node,
      const std::string& op_type,
      const std::string& input1,
      const std::string& input2,
      const std::string& output) {
    onnx::NodeProto* node = add_node();
    node->set_op_type(op_type);
    node->add_input(input1);
    node->add_input(input2);
    node->add_output(output);
    return node;
  }

  onnx::NodeProto* create_unary_node(
      AddNodeFn add_node,
      const std::string& op_type,
      const std::string& input,
      const std::string& output) {
    onnx::NodeProto* node = add_node();
    node->set_op_type(op_type);
    node->add_input(input);
    node->add_output(output);
    return node;
  }

  // GraphProto convenience wrappers
  onnx::NodeProto* create_binary_node(
      onnx::GraphProto* graph,
      const std::string& op_type,
      const std::string& input1,
      const std::string& input2,
      const std::string& output) {
    return create_binary_node([graph]() { return graph->add_node(); },
                               op_type, input1, input2, output);
  }

  onnx::NodeProto* create_unary_node(
      onnx::GraphProto* graph,
      const std::string& op_type,
      const std::string& input,
      const std::string& output) {
    return create_unary_node([graph]() { return graph->add_node(); },
                              op_type, input, output);
  }

  // Callback-based process_operation (main implementation)
  bool process_operation(
      AddNodeFn add_node,
      const Function& f,
      casadi_int op,
      casadi_int k,
      const std::vector<casadi_int>& i_vec,
      const std::vector<casadi_int>& o_vec,
      std::map<casadi_int, std::string>& work_to_onnx,
      const std::string& node_output) {

    onnx::NodeProto* node = nullptr;

    // Try simple operations using centralized lookup table
    const OpMapping* mapping = get_op_mapping(op);
    if (mapping) {
      if (mapping->arity == 1 && i_vec.size() >= 1 && o_vec.size() == 1) {
        create_unary_node(add_node, mapping->onnx_name, work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      } else if (mapping->arity == 2 && i_vec.size() >= 2 && o_vec.size() == 1) {
        create_binary_node(add_node, mapping->onnx_name, work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }
    }

    // Handle complex operations that can't use the simple lookup
    switch (op) {
      case OP_INPUT: {
        std::string input_name = f.name_in(i_vec[0]);
        if (input_name.empty()) {
          input_name = "input_" + std::to_string(i_vec[0]);
        }
        create_unary_node(add_node, "Identity", input_name, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_OUTPUT: {
        std::string output_name = f.name_out(o_vec[0]);
        if (output_name.empty()) {
          output_name = "output_" + std::to_string(o_vec[0]);
        }
        std::string input_onnx_name = work_to_onnx[i_vec[0]];
        create_unary_node(add_node, "Identity", input_onnx_name, output_name);

        // Note: Graph outputs are now added by add_graph_outputs() at the end
        // of function_to_graph() or at the end of load() for the main graph.
        // We only need to track the output name in work_to_onnx.
        work_to_onnx[o_vec[0]] = output_name;
        return true;
      }

      case OP_CONST: {
        MX mx_const = f.instruction_MX(k);
        DM dm_const = static_cast<DM>(mx_const);

        node = add_node();
        node->set_op_type("Constant");
        node->add_output(node_output);

        onnx::AttributeProto* attr = node->add_attribute();
        attr->set_name("value");
        onnx::TensorProto* tensor = attr->mutable_t();
        tensor->set_data_type(onnx::TensorProto::DOUBLE);

        tensor->add_dims(dm_const.size1());
        if (dm_const.size2() > 1) tensor->add_dims(dm_const.size2());

        for (casadi_int idx = 0; idx < dm_const.numel(); ++idx) {
          tensor->add_double_data(static_cast<double>(dm_const->at(idx)));
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_MTIMES: {
        // OP_MTIMES is a fused multiply-add: output = input[1] * input[2] + input[0]
        // This is a GEMM-like operation with 3 inputs
        std::string matmul_result = "matmul_" + std::to_string(k);
        create_binary_node(add_node, "MatMul", work_to_onnx[i_vec[1]], work_to_onnx[i_vec[2]], matmul_result);
        create_binary_node(add_node, "Add", matmul_result, work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_SQ: {
        // Square operation: output = input^2, implemented as Mul(input, input)
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[0]],
                          work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_TWICE: {
        // Double operation: output = 2*input
        // Create constant node with value 2
        std::string const_name = "const_2_" + std::to_string(k);
        onnx::NodeProto* const_node = add_node();
        const_node->set_op_type("Constant");
        const_node->add_output(const_name);
        onnx::AttributeProto* attr = const_node->add_attribute();
        attr->set_name("value");
        onnx::TensorProto* tensor = attr->mutable_t();
        tensor->set_data_type(onnx::TensorProto::DOUBLE);
        tensor->add_double_data(2.0);

        // Multiply input by 2
        create_binary_node(add_node, "Mul", const_name, work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_INV:
        // Reciprocal operation: output = 1/input
        create_unary_node(add_node, "Reciprocal", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_SOLVE: {
        // Linear solve operation: solve A*x = b for x
        // Extract transpose flag using info() pattern
        MX mx_solve = f.instruction_MX(k);
        Dict info = mx_solve.info();
        bool tr = info["tr"];  // Whether to solve A'*x = b instead

        // ONNX does not have a native linear solver operator
        // This would require matrix inversion or iterative solvers
        casadi_error("ONNX export: OP_SOLVE (linear solver) is not supported. "
                    "ONNX does not provide native linear algebra solvers. "
                    "Consider using explicit matrix operations or iterative methods. "
                    "Transpose flag was: " + std::string(tr ? "true" : "false"));
        return false;
      }

      // Reduction operations
      case OP_NORM1:
        create_unary_node(add_node, "ReduceL1", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_NORM2:
      case OP_NORMF:
        // Both L2 norm and Frobenius norm map to ReduceL2
        create_unary_node(add_node, "ReduceL2", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_NORMINF: {
        // Infinity norm: max(abs(x))
        std::string abs_result = "abs_" + std::to_string(k);
        create_unary_node(add_node, "Abs", work_to_onnx[i_vec[0]], abs_result);
        create_unary_node(add_node, "ReduceMax", abs_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_MMIN:
        create_unary_node(add_node, "ReduceMin", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_MMAX:
        create_unary_node(add_node, "ReduceMax", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      // Logical operations
      case OP_NOT:
        create_unary_node(add_node, "Not", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_OR:
        create_binary_node(add_node, "Or", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_AND:
        create_binary_node(add_node, "And", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_NE: {
        // Not equal: implemented as Not(Equal(...))
        std::string equal_result = "equal_" + std::to_string(k);
        create_binary_node(add_node, "Equal", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], equal_result);
        create_unary_node(add_node, "Not", equal_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_IF_ELSE_ZERO: {
        // Conditional: if (condition != 0) return value else return 0
        // Implemented as: Where(condition, value, 0)

        // Create constant zero
        std::string zero_name = "const_0_" + std::to_string(k);
        onnx::NodeProto* zero_node = add_node();
        zero_node->set_op_type("Constant");
        zero_node->add_output(zero_name);
        onnx::AttributeProto* attr = zero_node->add_attribute();
        attr->set_name("value");
        onnx::TensorProto* tensor = attr->mutable_t();
        tensor->set_data_type(onnx::TensorProto::DOUBLE);
        tensor->add_double_data(0.0);

        // Create Where node: Where(condition, value, zero)
        onnx::NodeProto* where_node = add_node();
        where_node->set_op_type("Where");
        where_node->add_input(work_to_onnx[i_vec[0]]);  // condition
        where_node->add_input(work_to_onnx[i_vec[1]]);  // value if true
        where_node->add_input(zero_name);               // value if false (zero)
        where_node->add_output(node_output);

        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Complex matrix operations
      case OP_DOT: {
        // Dot product: dot(a, b) = sum(a .* b)
        // Implemented as: ReduceSum(Mul(a, b))
        std::string mul_result = "mul_" + std::to_string(k);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], mul_result);
        create_unary_node(add_node, "ReduceSum", mul_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_BILIN: {
        // Bilinear form: input[1].' * input[0] * input[2]
        // Chain: Transpose(input[1]) -> MatMul with input[0] -> MatMul with input[2]
        std::string transpose_result = "transpose_" + std::to_string(k);
        std::string matmul1_result = "matmul1_" + std::to_string(k);

        create_unary_node(add_node, "Transpose", work_to_onnx[i_vec[1]], transpose_result);
        create_binary_node(add_node, "MatMul", transpose_result, work_to_onnx[i_vec[0]], matmul1_result);
        create_binary_node(add_node, "MatMul", matmul1_result, work_to_onnx[i_vec[2]], node_output);

        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_RANK1: {
        // Rank-1 update: input[0] + input[1] * input[2] * input[3].'
        // Chain: Transpose(input[3]) -> MatMul(input[2], result) -> MatMul(input[1], result) -> Add(input[0], result)
        std::string transpose_result = "transpose_" + std::to_string(k);
        std::string matmul1_result = "matmul1_" + std::to_string(k);
        std::string matmul2_result = "matmul2_" + std::to_string(k);

        create_unary_node(add_node, "Transpose", work_to_onnx[i_vec[3]], transpose_result);
        create_binary_node(add_node, "MatMul", work_to_onnx[i_vec[2]], transpose_result, matmul1_result);
        create_binary_node(add_node, "MatMul", work_to_onnx[i_vec[1]], matmul1_result, matmul2_result);
        create_binary_node(add_node, "Add", work_to_onnx[i_vec[0]], matmul2_result, node_output);

        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Tensor operations
      case OP_TRANSPOSE:
        create_unary_node(add_node, "Transpose", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_RESHAPE: {
        // Get the MX node to extract target shape info
        MX mx_reshape = f.instruction_MX(k);
        auto sp = mx_reshape.sparsity();

        node = add_node();
        node->set_op_type("Reshape");
        node->add_input(work_to_onnx[i_vec[0]]);

        // Create shape constant as second input
        std::string shape_name = "shape_" + std::to_string(k);
        onnx::NodeProto* shape_node = add_node();
        shape_node->set_op_type("Constant");
        shape_node->add_output(shape_name);

        onnx::AttributeProto* shape_attr = shape_node->add_attribute();
        shape_attr->set_name("value");
        onnx::TensorProto* shape_tensor = shape_attr->mutable_t();
        shape_tensor->set_data_type(onnx::TensorProto::INT64);
        shape_tensor->add_dims(2);
        shape_tensor->add_int64_data(sp.size1());
        shape_tensor->add_int64_data(sp.size2());

        node->add_input(shape_name);
        node->add_output(node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZCAT: {
        node = add_node();
        node->set_op_type("Concat");
        // Add all inputs
        for (casadi_int idx : i_vec) {
          node->add_input(work_to_onnx[idx]);
        }
        node->add_output(node_output);
        // Set axis attribute: axis=1 for horizontal concatenation (columns)
        onnx::AttributeProto* axis_attr = node->add_attribute();
        axis_attr->set_name("axis");
        axis_attr->set_i(1);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_VERTCAT: {
        node = add_node();
        node->set_op_type("Concat");
        // Add all inputs
        for (casadi_int idx : i_vec) {
          node->add_input(work_to_onnx[idx]);
        }
        node->add_output(node_output);
        // Set axis attribute: axis=0 for vertical concatenation (rows)
        onnx::AttributeProto* axis_attr = node->add_attribute();
        axis_attr->set_name("axis");
        axis_attr->set_i(0);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZSPLIT: {
        // Get the MX node and extract split info using info() pattern
        MX mx_split = f.instruction_MX(k);
        Dict info = mx_split.info();
        std::vector<casadi_int> offset = info["offset"];

        node = add_node();
        node->set_op_type("Split");
        node->add_input(work_to_onnx[i_vec[0]]);

        // Set axis attribute: axis=1 for horizontal split (columns)
        onnx::AttributeProto* axis_attr = node->add_attribute();
        axis_attr->set_name("axis");
        axis_attr->set_i(1);

        // Calculate split sizes from offset vector
        // offset = [0, size1, size1+size2, ...] -> split_sizes = [size1, size2, ...]
        std::vector<casadi_int> split_sizes;
        for (casadi_int j = 0; j < offset.size() - 1; ++j) {
          split_sizes.push_back(offset[j+1] - offset[j]);
        }

        // Add split sizes as attribute if they're not all equal
        bool all_equal = true;
        if (!split_sizes.empty()) {
          casadi_int first_size = split_sizes[0];
          for (casadi_int sz : split_sizes) {
            if (sz != first_size) {
              all_equal = false;
              break;
            }
          }
        }

        if (!all_equal) {
          onnx::AttributeProto* split_attr = node->add_attribute();
          split_attr->set_name("split");
          for (casadi_int sz : split_sizes) {
            split_attr->add_ints(sz);
          }
        }

        // Add all outputs
        for (casadi_int j = 0; j < o_vec.size(); ++j) {
          std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
          node->add_output(output_name);
          work_to_onnx[o_vec[j]] = output_name;
        }
        return true;
      }

      case OP_VERTSPLIT: {
        // Get the MX node and extract split info using info() pattern
        MX mx_split = f.instruction_MX(k);
        Dict info = mx_split.info();
        std::vector<casadi_int> offset = info["offset"];

        node = add_node();
        node->set_op_type("Split");
        node->add_input(work_to_onnx[i_vec[0]]);

        // Set axis attribute: axis=0 for vertical split (rows)
        onnx::AttributeProto* axis_attr = node->add_attribute();
        axis_attr->set_name("axis");
        axis_attr->set_i(0);

        // Calculate split sizes from offset vector
        // offset = [0, size1, size1+size2, ...] -> split_sizes = [size1, size2, ...]
        std::vector<casadi_int> split_sizes;
        for (casadi_int j = 0; j < offset.size() - 1; ++j) {
          split_sizes.push_back(offset[j+1] - offset[j]);
        }

        // Add split sizes as attribute if they're not all equal
        bool all_equal = true;
        if (!split_sizes.empty()) {
          casadi_int first_size = split_sizes[0];
          for (casadi_int sz : split_sizes) {
            if (sz != first_size) {
              all_equal = false;
              break;
            }
          }
        }

        if (!all_equal) {
          onnx::AttributeProto* split_attr = node->add_attribute();
          split_attr->set_name("split");
          for (casadi_int sz : split_sizes) {
            split_attr->add_ints(sz);
          }
        }

        // Add all outputs
        for (casadi_int j = 0; j < o_vec.size(); ++j) {
          std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
          node->add_output(output_name);
          work_to_onnx[o_vec[j]] = output_name;
        }
        return true;
      }

      // Operations not handled - caller should handle these
      case OP_CALL:
      case OP_SUBREF:
        return false;

      default:
        // Unknown operation - return false to let caller handle
        return false;
    }

    return false;
  }

  // GraphProto convenience wrapper for process_operation
  bool process_operation(
      onnx::GraphProto* graph,
      const Function& f,
      casadi_int op,
      casadi_int k,
      const std::vector<casadi_int>& i_vec,
      const std::vector<casadi_int>& o_vec,
      std::map<casadi_int, std::string>& work_to_onnx,
      const std::string& node_output) {
    return process_operation([graph]() { return graph->add_node(); },
                              f, op, k, i_vec, o_vec, work_to_onnx, node_output);
  }

} // namespace casadi
/// \endcond
