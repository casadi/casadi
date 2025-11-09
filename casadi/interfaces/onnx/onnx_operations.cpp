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


#include "onnx_model.hpp"

/// \cond INTERNAL
namespace casadi {

  // ========== Operation Mapping ==========
  //
  // Centralized mapping between CasADi opcodes and ONNX operation names.
  // Used by both export (CasADi → ONNX) and import (ONNX → CasADi).
  //
  // Single source of truth for all simple operations that have direct
  // CasADi ↔ ONNX equivalents.
  //
  // OpMapping struct is defined in onnx_model.hpp

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
    {OP_TRANSPOSE, "Transpose", 1},
    {OP_NORM1, "ReduceL1", 1},
    {OP_NORMF, "ReduceL2", 1},  // NORMF first: import uses norm_fro (works for vectors and matrices)
    {OP_NORM2, "ReduceL2", 1},
    {OP_MMIN, "ReduceMin", 1},
    {OP_MMAX, "ReduceMax", 1},
    // Binary operations
    {OP_ADD, "Add", 2},
    {OP_SUB, "Sub", 2},
    {OP_MUL, "Mul", 2},
    {OP_DIV, "Div", 2},
    {OP_POW, "Pow", 2},
    {OP_CONSTPOW, "Pow", 2},
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

  // Helper to create a Constant node with a tensor value
  // This ensures the attribute type is properly set (fixes ONNX Runtime compatibility)
  onnx::TensorProto* create_constant_tensor(
      AddNodeFn add_node,
      const std::string& output_name,
      onnx::TensorProto::DataType data_type) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Constant");
    node->add_output(output_name);
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name("value");
    attr->set_type(onnx::AttributeProto::TENSOR);  // Required for ONNX Runtime
    onnx::TensorProto* tensor = attr->mutable_t();
    tensor->set_data_type(data_type);
    return tensor;
  }

  // Add a Constant node holding an INT64 tensor (default shape: 1-D of data.size())
  void add_int_constant(AddNodeFn add_node, const std::string& name,
                        const std::vector<casadi_int>& data,
                        std::vector<casadi_int> dims = {}) {
    onnx::TensorProto* t = create_constant_tensor(add_node, name, onnx::TensorProto::INT64);
    if (dims.empty()) dims = {static_cast<casadi_int>(data.size())};
    for (casadi_int d : dims) t->add_dims(d);
    for (casadi_int v : data) t->add_int64_data(v);
  }

  // Add a Constant node holding a real tensor in the configured type (empty dims => scalar)
  void Onnx::add_real_constant(AddNodeFn add_node, const std::string& name,
                                         const std::vector<double>& data,
                                         const std::vector<casadi_int>& dims) {
    onnx::TensorProto* t = create_constant_tensor(add_node, name, real_type());
    for (casadi_int d : dims) t->add_dims(d);
    if (real_type() == onnx::TensorProto::FLOAT) {
      for (double v : data) t->add_float_data(static_cast<float>(v));
    } else {
      for (double v : data) t->add_double_data(v);
    }
  }

  // Add an integer attribute (e.g. axis) to a node
  void add_int_attribute(onnx::NodeProto* node, const std::string& name, casadi_int value) {
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name(name);
    attr->set_type(onnx::AttributeProto::INT);
    attr->set_i(value);
  }

  // Add an integer-list attribute (e.g. scan_input_axes) to a node
  void add_ints_attribute(onnx::NodeProto* node, const std::string& name,
                          const std::vector<casadi_int>& values) {
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name(name);
    attr->set_type(onnx::AttributeProto::INTS);
    for (casadi_int v : values) attr->add_ints(v);
  }

  // GraphProto convenience wrapper for add_int_constant
  void add_int_constant(onnx::GraphProto* graph, const std::string& name,
                        const std::vector<casadi_int>& data, std::vector<casadi_int> dims) {
    add_int_constant([graph]() { return graph->add_node(); }, name, data, dims);
  }

  // Gather(data, indices) along an axis -> output
  onnx::NodeProto* create_gather_node(AddNodeFn add_node, const std::string& data,
                                      const std::string& indices, const std::string& output,
                                      casadi_int axis = 0) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Gather");
    node->add_input(data);
    node->add_input(indices);
    node->add_output(output);
    add_int_attribute(node, "axis", axis);
    return node;
  }

  // Slice(data) on the given axes from starts to ends (step 1), emitting the index constants
  onnx::NodeProto* create_slice_node(AddNodeFn add_node, const std::string& uniq,
      const std::string& data, const std::vector<casadi_int>& starts,
      const std::vector<casadi_int>& ends, const std::vector<casadi_int>& axes,
      const std::string& output) {
    add_int_constant(add_node, uniq + "_starts", starts);
    add_int_constant(add_node, uniq + "_ends", ends);
    add_int_constant(add_node, uniq + "_axes", axes);
    onnx::NodeProto* node = add_node();
    node->set_op_type("Slice");
    node->add_input(data);
    node->add_input(uniq + "_starts");
    node->add_input(uniq + "_ends");
    node->add_input(uniq + "_axes");
    node->add_output(output);
    return node;
  }

  // Reshape `data` to `dims` with CasADi column-major semantics. ONNX Reshape is row-major,
  // so for a genuine 2-D target we emulate column-major as Transpose . Reshape(reverse) .
  // Transpose; for <=1-D targets column- and row-major coincide, so a plain Reshape suffices.
  void emit_colmajor_reshape(AddNodeFn add_node, const std::string& data,
                             const std::vector<casadi_int>& dims, const std::string& output,
                             const std::string& uniq) {
    if (dims.size() < 2) {
      add_int_constant(add_node, uniq + "_s", dims);
      create_unary_node(add_node, "Reshape", data, output)->add_input(uniq + "_s");
    } else {
      std::string t1 = uniq + "_t1";
      create_unary_node(add_node, "Transpose", data, t1);
      add_int_constant(add_node, uniq + "_s", {dims[1], dims[0]});
      std::string r = uniq + "_r";
      create_unary_node(add_node, "Reshape", t1, r)->add_input(uniq + "_s");
      create_unary_node(add_node, "Transpose", r, output);
    }
  }

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

  // Cast(input) -> output of the given ONNX data type
  onnx::NodeProto* create_cast_node(
      AddNodeFn add_node,
      const std::string& input,
      const std::string& output,
      onnx::TensorProto::DataType to_type) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Cast");
    node->add_input(input);
    node->add_output(output);
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name("to");
    attr->set_type(onnx::AttributeProto::INT);
    attr->set_i(to_type);
    return node;
  }

  // Gemm: output = (A or A')*(B or B') [+ C]. C empty -> 2-input form (alpha=beta=1)
  onnx::NodeProto* create_gemm_node(
      AddNodeFn add_node,
      const std::string& A, const std::string& B, const std::string& C,
      const std::string& output, bool transA = false, bool transB = false) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Gemm");
    node->add_input(A);
    node->add_input(B);
    if (!C.empty()) node->add_input(C);
    node->add_output(output);
    if (transA) add_int_attribute(node, "transA", 1);
    if (transB) add_int_attribute(node, "transB", 1);
    return node;
  }

  // Where(condition, if_true, if_false) -> output
  onnx::NodeProto* create_where_node(
      AddNodeFn add_node,
      const std::string& cond,
      const std::string& if_true,
      const std::string& if_false,
      const std::string& output) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Where");
    node->add_input(cond);
    node->add_input(if_true);
    node->add_input(if_false);
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
  bool Onnx::process_operation(
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
        std::string input_name = onnx_input_name(f, i_vec[0]);

        // Get MX info to check if we're accessing a specific element of a larger input
        MX mx_input = f.instruction_MX(k);
        Dict info = mx_input.info();
        casadi_int offset = info["offset"];
        casadi_int input_numel = f.numel_in(i_vec[0]);
        casadi_int output_numel = mx_input.numel();

        if (output_numel < input_numel) {
          // Multi-element input: use Gather to extract the specific element at offset
          std::string index_name = "input_idx_" + std::to_string(k);
          add_int_constant(add_node, index_name, {offset});

          // Create Gather node to extract element at offset (axis=0)
          node = add_node();
          node->set_op_type("Gather");
          node->add_input(input_name);
          node->add_input(index_name);
          node->add_output(node_output);
          add_int_attribute(node, "axis", 0);
        } else {
          // Full input access (output size matches input size), use Identity
          create_unary_node(add_node, "Identity", input_name, node_output);
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_OUTPUT: {
        // Connect the work value to the output name; o_vec[0] is the output slot, not a work index
        create_unary_node(add_node, "Identity", work_to_onnx[i_vec[0]], onnx_output_name(f, o_vec[0]));
        return true;
      }

      case OP_CONST: {
        // Always emit 2-D dims so elementwise broadcasting matches CasADi's 2-D convention
        // (a rank-1 [n] would broadcast (n,1)+(n,) -> (n,n)).
        DM dm_const = static_cast<DM>(f.instruction_MX(k));
        std::vector<double> data(dm_const->begin(), dm_const->end());
        add_real_constant(add_node, node_output, data, {dm_const.size1(), dm_const.size2()});
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_MTIMES: {
        // Fused multiply-add input[1]*input[2] + input[0] -> a single Gemm
        create_gemm_node(add_node, work_to_onnx[i_vec[1]], work_to_onnx[i_vec[2]],
                         work_to_onnx[i_vec[0]], node_output);
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
        std::string const_name = "const_2_" + std::to_string(k);
        add_real_constant(add_node, const_name, {2.0});
        create_binary_node(add_node, "Mul", const_name, work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_DETERMINANT:
        create_unary_node(add_node, "Det", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_LOGSUMEXP:
        create_unary_node(add_node, "ReduceLogSumExp", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_LOG1P: {
        // log(1 + x)
        std::string one = "one_" + std::to_string(k), s = "log1p_" + std::to_string(k);
        add_real_constant(add_node, one, {1.0});
        create_binary_node(add_node, "Add", one, work_to_onnx[i_vec[0]], s);
        create_unary_node(add_node, "Log", s, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_EXPM1: {
        // exp(x) - 1
        std::string e = "exp_" + std::to_string(k), one = "one_" + std::to_string(k);
        create_unary_node(add_node, "Exp", work_to_onnx[i_vec[0]], e);
        add_real_constant(add_node, one, {1.0});
        create_binary_node(add_node, "Sub", e, one, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HYPOT: {
        // sqrt(x^2 + y^2)
        std::string x2 = "hx_" + std::to_string(k), y2 = "hy_" + std::to_string(k);
        std::string s = "hs_" + std::to_string(k);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[0]], x2);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[1]], work_to_onnx[i_vec[1]], y2);
        create_binary_node(add_node, "Add", x2, y2, s);
        create_unary_node(add_node, "Sqrt", s, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Pure pass-throughs (drop CasADi-internal markers)
      case OP_ASSIGN:
      case OP_LIFT:
        create_unary_node(add_node, "Identity", work_to_onnx[i_vec[0]], node_output);
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

      case OP_NORMINF: {
        // Infinity norm: max(abs(x))
        std::string abs_result = "abs_" + std::to_string(k);
        create_unary_node(add_node, "Abs", work_to_onnx[i_vec[0]], abs_result);
        create_unary_node(add_node, "ReduceMax", abs_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Comparisons output a bool tensor in ONNX; cast to double for the double-typed graph
      case OP_LT: {
        std::string b = "cmp_" + std::to_string(k);
        create_binary_node(add_node, "Less", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], b);
        create_cast_node(add_node, b, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_LE: {
        std::string b = "cmp_" + std::to_string(k);
        create_binary_node(add_node, "LessOrEqual", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], b);
        create_cast_node(add_node, b, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_EQ: {
        std::string b = "cmp_" + std::to_string(k);
        create_binary_node(add_node, "Equal", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], b);
        create_cast_node(add_node, b, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Logic ops need bool operands: cast 0/1 doubles to bool, apply, cast back
      case OP_NOT: {
        std::string b = "bin_" + std::to_string(k);
        create_cast_node(add_node, work_to_onnx[i_vec[0]], b, onnx::TensorProto::BOOL);
        std::string r = "bres_" + std::to_string(k);
        create_unary_node(add_node, "Not", b, r);
        create_cast_node(add_node, r, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_AND:
      case OP_OR: {
        std::string b0 = "bin0_" + std::to_string(k), b1 = "bin1_" + std::to_string(k);
        create_cast_node(add_node, work_to_onnx[i_vec[0]], b0, onnx::TensorProto::BOOL);
        create_cast_node(add_node, work_to_onnx[i_vec[1]], b1, onnx::TensorProto::BOOL);
        std::string r = "bres_" + std::to_string(k);
        create_binary_node(add_node, op == OP_AND ? "And" : "Or", b0, b1, r);
        create_cast_node(add_node, r, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_FMIN:
        create_binary_node(add_node, "Min", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_FMOD: {
        node = add_node();
        node->set_op_type("Mod");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(work_to_onnx[i_vec[1]]);
        node->add_output(node_output);
        add_int_attribute(node, "fmod", 1);  // fmod=1 required for floating point types
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_COPYSIGN: {
        // copysign(x, y) = sign(y) * abs(x)
        std::string sign_result = "sign_" + std::to_string(k);
        std::string abs_result = "abs_" + std::to_string(k);
        create_unary_node(add_node, "Sign", work_to_onnx[i_vec[1]], sign_result);
        create_unary_node(add_node, "Abs", work_to_onnx[i_vec[0]], abs_result);
        create_binary_node(add_node, "Mul", sign_result, abs_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_FMAX:
        create_binary_node(add_node, "Max", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_NE: {
        // Not equal: Not(Equal(...)), then cast the bool result to double
        std::string equal_result = "equal_" + std::to_string(k);
        create_binary_node(add_node, "Equal", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], equal_result);
        std::string r = "bres_" + std::to_string(k);
        create_unary_node(add_node, "Not", equal_result, r);
        create_cast_node(add_node, r, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_IF_ELSE_ZERO: {
        // if (condition != 0) value else 0  =>  Where(bool(condition), value, 0)
        std::string cond = "cond_" + std::to_string(k);
        create_cast_node(add_node, work_to_onnx[i_vec[0]], cond, onnx::TensorProto::BOOL);
        std::string zero_name = "const_0_" + std::to_string(k);
        add_real_constant(add_node, zero_name, {0.0});
        create_where_node(add_node, cond, work_to_onnx[i_vec[1]], zero_name, node_output);
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
        // Bilinear form input[1]' * input[0] * input[2]. Gemm(transA) avoids a separate
        // Transpose (which ORT would fuse into a float64-less FusedMatMul).
        std::string xa = "bilin_" + std::to_string(k);
        create_gemm_node(add_node, work_to_onnx[i_vec[1]], work_to_onnx[i_vec[0]], "", xa, true);
        create_binary_node(add_node, "MatMul", xa, work_to_onnx[i_vec[2]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_RANK1: {
        // Rank-1 update input[0] + input[1]*input[2]*input[3]'. Gemm(transB) forms x*y'
        // without a separate Transpose; input[1] is a scalar alpha.
        std::string xyt = "rank1_" + std::to_string(k);
        create_gemm_node(add_node, work_to_onnx[i_vec[2]], work_to_onnx[i_vec[3]], "", xyt, false, true);
        std::string scaled = "rank1s_" + std::to_string(k);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[1]], xyt, scaled);
        create_binary_node(add_node, "Add", work_to_onnx[i_vec[0]], scaled, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_EINSTEIN: {
        // Tensor contraction C_c += A_a * B_b. deps [C, A, B], operands column-vec flattened.
        // Reshape each operand (column-major) to its logical dims and contract with Einsum
        // using the labels in natural order; reshape the result back to a column vector.
        MX mx_e = f.instruction_MX(k);
        Dict info = mx_e.info();
        std::vector<casadi_int> la = info["a"], lb = info["b"], lc = info["c"];
        std::vector<casadi_int> da = info["dim_a"], db = info["dim_b"], dc = info["dim_c"];
        casadi_assert(da.size() <= 2 && db.size() <= 2 && dc.size() <= 2,
          "ONNX export: einstein with >2-index operands is not supported.");

        // Assign a letter to each distinct label (natural order)
        std::map<casadi_int, char> letter;
        char nxt = 'a';
        for (const std::vector<casadi_int>& labs : {la, lb, lc}) {
          for (casadi_int L : labs) if (!letter.count(L)) letter[L] = nxt++;
        }
        std::string sa, sb, sc;
        for (casadi_int L : la) sa += letter[L];
        for (casadi_int L : lb) sb += letter[L];
        for (casadi_int L : lc) sc += letter[L];

        std::string a_re = "ein_a_" + std::to_string(k), b_re = "ein_b_" + std::to_string(k);
        emit_colmajor_reshape(add_node, work_to_onnx[i_vec[1]], da, a_re, a_re);
        emit_colmajor_reshape(add_node, work_to_onnx[i_vec[2]], db, b_re, b_re);

        std::string ein_out = "ein_o_" + std::to_string(k);
        node = add_node();
        node->set_op_type("Einsum");
        node->add_input(a_re);
        node->add_input(b_re);
        node->add_output(ein_out);
        onnx::AttributeProto* eq = node->add_attribute();
        eq->set_name("equation");
        eq->set_type(onnx::AttributeProto::STRING);
        eq->set_s(sa + "," + sb + "->" + sc);

        // Flatten the contraction result back to a column vector and add the C accumulator
        casadi_int prod_c = 1;
        for (casadi_int d : dc) prod_c *= d;
        std::string ein_vec = "ein_v_" + std::to_string(k);
        emit_colmajor_reshape(add_node, ein_out, {prod_c, 1}, ein_vec, ein_vec);
        create_binary_node(add_node, "Add", work_to_onnx[i_vec[0]], ein_vec, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Tensor operations
      case OP_RESHAPE: {
        auto sp = f.instruction_MX(k).sparsity();
        emit_colmajor_reshape(add_node, work_to_onnx[i_vec[0]], {sp.size1(), sp.size2()},
                              node_output, "rs_" + std::to_string(k));
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZCAT:
      case OP_VERTCAT: {
        // axis=1 (columns) for horzcat, axis=0 (rows) for vertcat
        node = add_node();
        node->set_op_type("Concat");
        for (casadi_int idx : i_vec) node->add_input(work_to_onnx[idx]);
        node->add_output(node_output);
        add_int_attribute(node, "axis", op == OP_HORZCAT ? 1 : 0);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZSPLIT:
      case OP_VERTSPLIT: {
        // horzsplit -> axis=1, divide NNZ offsets by nrows; vertsplit -> axis=0, by ncols
        MX mx_split = f.instruction_MX(k);
        std::vector<casadi_int> offset = mx_split.info()["offset"];
        bool horz = (op == OP_HORZSPLIT);
        casadi_int stride = horz ? mx_split.dep(0).size1() : mx_split.dep(0).size2();

        std::vector<casadi_int> split_sizes;
        for (casadi_int j = 0; j + 1 < offset.size(); ++j) {
          split_sizes.push_back((offset[j+1] - offset[j]) / stride);
        }

        std::string split_sizes_name = "split_sizes_" + std::to_string(k);
        add_int_constant(add_node, split_sizes_name, split_sizes);

        node = add_node();
        node->set_op_type("Split");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(split_sizes_name);
        add_int_attribute(node, "axis", horz ? 1 : 0);

        for (casadi_int j = 0; j < o_vec.size(); ++j) {
          std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
          node->add_output(output_name);
          work_to_onnx[o_vec[j]] = output_name;
        }
        return true;
      }

      case OP_HORZREPMAT: {
        // repmat(x, 1, n) -> ONNX Tile with repeats [1, n]
        MX mx_repmat = f.instruction_MX(k);
        casadi_int n = mx_repmat.size2() / mx_repmat.dep(0).size2();

        std::string repeats_name = "repeats_" + std::to_string(k);
        add_int_constant(add_node, repeats_name, {1, n});

        node = add_node();
        node->set_op_type("Tile");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(repeats_name);
        node->add_output(node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZREPSUM: {
        // repsum(x, 1, n): inverse of repmat - split [rows, cols*n] into n parts and sum
        MX mx_repsum = f.instruction_MX(k);
        casadi_int output_cols = mx_repsum.size2();
        casadi_int n = mx_repsum.dep(0).size2() / output_cols;

        std::string split_sizes_name = "split_sizes_" + std::to_string(k);
        add_int_constant(add_node, split_sizes_name, std::vector<casadi_int>(n, output_cols));

        node = add_node();
        node->set_op_type("Split");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(split_sizes_name);
        add_int_attribute(node, "axis", 1);

        // Sum all n parts in a single variadic Sum node
        onnx::NodeProto* sum_node = add_node();
        sum_node->set_op_type("Sum");
        for (casadi_int i = 0; i < n; ++i) {
          std::string out_name = "split_" + std::to_string(k) + "_" + std::to_string(i);
          node->add_output(out_name);
          sum_node->add_input(out_name);
        }
        sum_node->add_output(node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_GETNONZEROS: {
        // x.nz[indices] in CasADi column-major order. ONNX indexing is row-major, so we
        // column-major-flatten the input, Gather/Slice the flat (a column vector), then
        // reshape the result back column-major. The 2-D 'inner'/'outer' form encodes flat
        // indices too (inner = first block's flat range, outer = block offsets).
        MX mx_getnonzeros = f.instruction_MX(k);
        Dict info = mx_getnonzeros.info();
        std::string data = work_to_onnx[i_vec[0]];
        std::string uniq = std::to_string(k);

        // Collect the column-major flat indices (or a contiguous unit-step range)
        std::vector<casadi_int> idx;
        bool contiguous = false;
        casadi_int c_start = 0, c_stop = 0;
        if (info.count("nz")) {
          idx = info["nz"];
        } else if (info.count("slice")) {
          Dict s = info["slice"];
          casadi_int start = s["start"], stop = s["stop"], step = s["step"];
          if (step == 1) { contiguous = true; c_start = start; c_stop = stop; }
          else for (casadi_int i = start; i < stop; i += step) idx.push_back(i);
        } else if (info.count("inner")) {
          Dict inner_info = info["inner"], outer_info = info["outer"];
          casadi_int is = inner_info["start"], ip = inner_info["stop"], ist = inner_info["step"];
          casadi_int os = outer_info["start"], op2 = outer_info["stop"], ost = outer_info["step"];
          for (casadi_int o = os; o < op2; o += ost) {
            for (casadi_int i = is; i < ip; i += ist) idx.push_back(o + i);
          }
        } else {
          return false;  // unknown variant
        }

        // Flatten the input column-major (skip when it is already a column vector)
        std::string flat = data;
        if (mx_getnonzeros.dep(0).size2() != 1) {
          flat = "gnz_flat_" + uniq;
          emit_colmajor_reshape(add_node, data, {mx_getnonzeros.dep(0).numel(), 1}, flat,
                                "gnz_flat_" + uniq);
        }

        // Gather/Slice the flat; write straight to the output if it is the (k,1) column
        bool out_is_col = (mx_getnonzeros.size2() == 1);
        std::string gathered = out_is_col ? node_output : ("gnz_g_" + uniq);
        if (contiguous) {
          create_slice_node(add_node, uniq, flat, {c_start}, {c_stop}, {0}, gathered);
        } else {
          add_int_constant(add_node, "indices_" + uniq, idx);
          create_gather_node(add_node, flat, "indices_" + uniq, gathered);
        }
        if (!out_is_col) {
          emit_colmajor_reshape(add_node, gathered,
                                {mx_getnonzeros.size1(), mx_getnonzeros.size2()}, node_output,
                                "gnz_out_" + uniq);
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_SETNONZEROS: {
        // y = x with y.nz[idx] = z (column-major). Flatten x column-major, ScatterElements on
        // the flat (axis 0), reshape back -- mirrors GETNONZEROS.
        MX mx_set = f.instruction_MX(k);
        Dict info = mx_set.info();
        std::string uniq = std::to_string(k);
        bool add = info["add"];
        casadi_assert(!add, "ONNX export: scatter-add (setnonzeros with add) is not supported.");

        std::vector<casadi_int> idx;
        if (info.count("nz")) {
          idx = info["nz"];
        } else if (info.count("slice")) {
          Dict s = info["slice"];
          casadi_int start = s["start"], stop = s["stop"], step = s["step"];
          for (casadi_int i = start; i < stop; i += step) idx.push_back(i);
        } else if (info.count("inner")) {
          Dict ii = info["inner"], oo = info["outer"];
          casadi_int is = ii["start"], ip = ii["stop"], ist = ii["step"];
          casadi_int os = oo["start"], op2 = oo["stop"], ost = oo["step"];
          for (casadi_int o = os; o < op2; o += ost) {
            for (casadi_int i = is; i < ip; i += ist) idx.push_back(o + i);
          }
        } else {
          casadi_error("ONNX export: unsupported setnonzeros pattern.");
        }

        // Flatten data and updates column-major (skip when already a column vector)
        bool data_col = (mx_set.dep(0).size2() == 1);
        std::string flat = work_to_onnx[i_vec[0]];
        if (!data_col) {
          flat = "snz_flat_" + uniq;
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[0]], {mx_set.dep(0).numel(), 1},
                                flat, "snz_flat_" + uniq);
        }
        std::string upd = work_to_onnx[i_vec[1]];
        if (mx_set.dep(1).size2() != 1) {
          upd = "snz_upd_" + uniq;
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[1]], {mx_set.dep(1).numel(), 1},
                                upd, "snz_upd_" + uniq);
        }

        std::string ind = "snz_idx_" + uniq;
        add_int_constant(add_node, ind, idx, {static_cast<casadi_int>(idx.size()), 1});
        std::string scattered = data_col ? node_output : ("snz_s_" + uniq);
        node = add_node();
        node->set_op_type("ScatterElements");
        node->add_input(flat);
        node->add_input(ind);
        node->add_input(upd);
        node->add_output(scattered);
        add_int_attribute(node, "axis", 0);

        if (!data_col) {
          emit_colmajor_reshape(add_node, scattered, {mx_set.dep(0).size1(), mx_set.dep(0).size2()},
                                node_output, "snz_out_" + uniq);
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_ATAN2: {
        // atan2(y, x) = atan(y/x) + correction for x<0, with origin→0
        std::string y = work_to_onnx[i_vec[0]], x = work_to_onnx[i_vec[1]];
        std::string p = "atan2_" + std::to_string(k) + "_";

        add_real_constant(add_node, p+"c0", {0.0});
        add_real_constant(add_node, p+"pi", {M_PI});
        create_binary_node(add_node, "Div", y, x, p+"r");
        create_unary_node(add_node, "Atan", p+"r", p+"b");
        create_binary_node(add_node, "Less", x, p+"c0", p+"xn");
        create_binary_node(add_node, "Less", y, p+"c0", p+"yn");
        create_unary_node(add_node, "Neg", p+"pi", p+"npi");
        create_where_node(add_node, p+"yn", p+"npi", p+"pi", p+"corr"); // x<0: -π if y<0 else +π
        create_where_node(add_node, p+"xn", p+"corr", p+"c0", p+"adj"); // 0 when x>=0
        create_binary_node(add_node, "Add", p+"b", p+"adj", p+"res");
        create_binary_node(add_node, "Equal", x, p+"c0", p+"xz");
        create_binary_node(add_node, "Equal", y, p+"c0", p+"yz");
        create_binary_node(add_node, "And", p+"xz", p+"yz", p+"oz");
        create_where_node(add_node, p+"oz", p+"c0", p+"res", node_output); // origin → 0

        work_to_onnx[o_vec[0]] = node_output;
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
  bool Onnx::process_operation(
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
