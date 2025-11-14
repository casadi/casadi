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

    // Binary operations
    if (op_type == "Add") {
      casadi_assert(node_inputs.size() >= 2, "Add requires 2 inputs");
      output = node_inputs[0] + node_inputs[1];
    } else if (op_type == "Sub") {
      casadi_assert(node_inputs.size() >= 2, "Sub requires 2 inputs");
      output = node_inputs[0] - node_inputs[1];
    } else if (op_type == "Mul") {
      casadi_assert(node_inputs.size() >= 2, "Mul requires 2 inputs");
      output = node_inputs[0] * node_inputs[1];
    } else if (op_type == "Div") {
      casadi_assert(node_inputs.size() >= 2, "Div requires 2 inputs");
      output = node_inputs[0] / node_inputs[1];
    } else if (op_type == "Pow") {
      casadi_assert(node_inputs.size() >= 2, "Pow requires 2 inputs");
      output = pow(node_inputs[0], node_inputs[1]);

    // Trigonometric
    } else if (op_type == "Sin") {
      casadi_assert(node_inputs.size() >= 1, "Sin requires 1 input");
      output = sin(node_inputs[0]);
    } else if (op_type == "Cos") {
      casadi_assert(node_inputs.size() >= 1, "Cos requires 1 input");
      output = cos(node_inputs[0]);
    } else if (op_type == "Tan") {
      casadi_assert(node_inputs.size() >= 1, "Tan requires 1 input");
      output = tan(node_inputs[0]);

    // Inverse Trigonometric
    } else if (op_type == "Asin") {
      casadi_assert(node_inputs.size() >= 1, "Asin requires 1 input");
      output = asin(node_inputs[0]);
    } else if (op_type == "Acos") {
      casadi_assert(node_inputs.size() >= 1, "Acos requires 1 input");
      output = acos(node_inputs[0]);
    } else if (op_type == "Atan") {
      casadi_assert(node_inputs.size() >= 1, "Atan requires 1 input");
      output = atan(node_inputs[0]);

    // Hyperbolic
    } else if (op_type == "Sinh") {
      casadi_assert(node_inputs.size() >= 1, "Sinh requires 1 input");
      output = sinh(node_inputs[0]);
    } else if (op_type == "Cosh") {
      casadi_assert(node_inputs.size() >= 1, "Cosh requires 1 input");
      output = cosh(node_inputs[0]);
    } else if (op_type == "Tanh") {
      casadi_assert(node_inputs.size() >= 1, "Tanh requires 1 input");
      output = tanh(node_inputs[0]);

    // Inverse Hyperbolic
    } else if (op_type == "Asinh") {
      casadi_assert(node_inputs.size() >= 1, "Asinh requires 1 input");
      output = asinh(node_inputs[0]);
    } else if (op_type == "Acosh") {
      casadi_assert(node_inputs.size() >= 1, "Acosh requires 1 input");
      output = acosh(node_inputs[0]);
    } else if (op_type == "Atanh") {
      casadi_assert(node_inputs.size() >= 1, "Atanh requires 1 input");
      output = atanh(node_inputs[0]);

    // Exponential/Logarithmic
    } else if (op_type == "Exp") {
      casadi_assert(node_inputs.size() >= 1, "Exp requires 1 input");
      output = exp(node_inputs[0]);
    } else if (op_type == "Log") {
      casadi_assert(node_inputs.size() >= 1, "Log requires 1 input");
      output = log(node_inputs[0]);
    } else if (op_type == "Sqrt") {
      casadi_assert(node_inputs.size() >= 1, "Sqrt requires 1 input");
      output = sqrt(node_inputs[0]);

    // Rounding
    } else if (op_type == "Ceil") {
      casadi_assert(node_inputs.size() >= 1, "Ceil requires 1 input");
      output = ceil(node_inputs[0]);
    } else if (op_type == "Floor") {
      casadi_assert(node_inputs.size() >= 1, "Floor requires 1 input");
      output = floor(node_inputs[0]);

    // Other
    } else if (op_type == "Abs") {
      casadi_assert(node_inputs.size() >= 1, "Abs requires 1 input");
      output = fabs(node_inputs[0]);
    } else if (op_type == "Sign") {
      casadi_assert(node_inputs.size() >= 1, "Sign requires 1 input");
      output = sign(node_inputs[0]);
    } else if (op_type == "Neg") {
      casadi_assert(node_inputs.size() >= 1, "Neg requires 1 input");
      output = -node_inputs[0];
    } else if (op_type == "Erf") {
      casadi_assert(node_inputs.size() >= 1, "Erf requires 1 input");
      output = erf(node_inputs[0]);

    // Utilities
    } else if (op_type == "Identity") {
      casadi_assert(node_inputs.size() >= 1, "Identity requires 1 input");
      output = node_inputs[0];

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

    } else {
      casadi_error("Unsupported operation '" + op_type + "'");
    }

    return output;
  }

} // namespace casadi
/// \endcond
