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

  void OnnxTranslator::load(const Function& f) {
    // Create ONNX model from CasADi Function
    model_.Clear();
    model_.set_ir_version(8);

    // Set producer info
    model_.set_producer_name("CasADi");
    model_.set_producer_version(CASADI_VERSION_STRING);

    // Add opset import (required by ONNX spec)
    onnx::OperatorSetIdProto* opset = model_.add_opset_import();
    opset->set_domain("");  // Empty string means default ONNX domain
    opset->set_version(13);  // Use opset version 13

    // Create graph
    onnx::GraphProto* graph = model_.mutable_graph();
    graph->set_name(f.name());

    // Add inputs to graph
    for (casadi_int i = 0; i < f.n_in(); ++i) {
      onnx::ValueInfoProto* input = graph->add_input();
      std::string input_name = f.name_in(i);
      if (input_name.empty()) {
        input_name = "input_" + std::to_string(i);
      }
      input->set_name(input_name);

      // Set tensor type and shape
      onnx::TypeProto* type = input->mutable_type();
      onnx::TypeProto::Tensor* tensor_type = type->mutable_tensor_type();
      tensor_type->set_elem_type(onnx::TensorProto::DOUBLE);

      // Add shape dimensions
      onnx::TensorShapeProto* shape = tensor_type->mutable_shape();
      auto sp = f.sparsity_in(i);
      shape->add_dim()->set_dim_value(sp.size1());
      shape->add_dim()->set_dim_value(sp.size2());
    }

    // Walk through the algorithm and build ONNX nodes
    // Map from work vector index to the ONNX node name that produced it
    std::map<casadi_int, std::string> work_to_onnx;

    casadi_int n_instr = f.n_instructions();

    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i = f.instruction_input(k);

      // Node name for this operation's result - use instruction index for uniqueness
      std::string node_output = "n" + std::to_string(k);

      // Create ONNX node based on operation type
      onnx::NodeProto* node = nullptr;

      switch (op) {
        case OP_INPUT:
          // Input nodes are already added, map to input names
          {
            std::string input_name = f.name_in(i[0]);
            if (input_name.empty()) {
              input_name = "input_" + std::to_string(i[0]);
            }
            // Create Identity node to map input to work vector
            node = graph->add_node();
            node->set_op_type("Identity");
            node->add_input(input_name);
            node->add_output(node_output);
            // Track this output in the work vector mapping
            work_to_onnx[o[0]] = node_output;
          }
          break;

        case OP_OUTPUT:
          // Output nodes: map from work vector to output names
          {
            std::string output_name = f.name_out(o[0]);
            if (output_name.empty()) {
              output_name = "output_" + std::to_string(o[0]);
            }
            // Get the ONNX name that produced this work vector value
            std::string input_onnx_name = work_to_onnx[i[0]];

            // Create Identity node to map work vector to output
            node = graph->add_node();
            node->set_op_type("Identity");
            node->add_input(input_onnx_name);
            node->add_output(output_name);

            // Add this as a graph output
            onnx::ValueInfoProto* graph_output = graph->add_output();
            graph_output->set_name(output_name);
            onnx::TypeProto* type = graph_output->mutable_type();
            onnx::TypeProto::Tensor* tensor_type = type->mutable_tensor_type();
            tensor_type->set_elem_type(onnx::TensorProto::DOUBLE);
            // Add shape from function output sparsity
            onnx::TensorShapeProto* shape = tensor_type->mutable_shape();
            auto sp = f.sparsity_out(o[0]);
            shape->add_dim()->set_dim_value(sp.size1());
            shape->add_dim()->set_dim_value(sp.size2());
          }
          break;

        case OP_CONST: {
          // Constant node - convert MX to DM
          MX mx_const = f.instruction_MX(k);
          DM dm_const = static_cast<DM>(mx_const);

          node = graph->add_node();
          node->set_op_type("Constant");
          node->add_output(node_output);

          onnx::AttributeProto* attr = node->add_attribute();
          attr->set_name("value");
          onnx::TensorProto* tensor = attr->mutable_t();
          tensor->set_data_type(onnx::TensorProto::DOUBLE);

          // Add shape
          tensor->add_dims(dm_const.size1());
          if (dm_const.size2() > 1) tensor->add_dims(dm_const.size2());

          // Add data
          for (casadi_int idx = 0; idx < dm_const.numel(); ++idx) {
            tensor->add_double_data(static_cast<double>(dm_const->at(idx)));
          }
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_ADD: {
          node = graph->add_node();
          node->set_op_type("Add");
          node->add_input(work_to_onnx[i[0]]);
          node->add_input(work_to_onnx[i[1]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_SUB: {
          node = graph->add_node();
          node->set_op_type("Sub");
          node->add_input(work_to_onnx[i[0]]);
          node->add_input(work_to_onnx[i[1]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_MUL: {
          node = graph->add_node();
          node->set_op_type("Mul");
          node->add_input(work_to_onnx[i[0]]);
          node->add_input(work_to_onnx[i[1]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_DIV: {
          node = graph->add_node();
          node->set_op_type("Div");
          node->add_input(work_to_onnx[i[0]]);
          node->add_input(work_to_onnx[i[1]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_SIN: {
          node = graph->add_node();
          node->set_op_type("Sin");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_COS: {
          node = graph->add_node();
          node->set_op_type("Cos");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_TAN: {
          node = graph->add_node();
          node->set_op_type("Tan");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_EXP: {
          node = graph->add_node();
          node->set_op_type("Exp");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_LOG: {
          node = graph->add_node();
          node->set_op_type("Log");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_SQRT: {
          node = graph->add_node();
          node->set_op_type("Sqrt");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_NEG: {
          node = graph->add_node();
          node->set_op_type("Neg");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_TANH: {
          node = graph->add_node();
          node->set_op_type("Tanh");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        // ========== Inverse Trigonometric ==========
        case OP_ASIN: {
          node = graph->add_node();
          node->set_op_type("Asin");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_ACOS: {
          node = graph->add_node();
          node->set_op_type("Acos");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_ATAN: {
          node = graph->add_node();
          node->set_op_type("Atan");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        // ========== Hyperbolic ==========
        case OP_SINH: {
          node = graph->add_node();
          node->set_op_type("Sinh");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_COSH: {
          node = graph->add_node();
          node->set_op_type("Cosh");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        // ========== Inverse Hyperbolic ==========
        case OP_ASINH: {
          node = graph->add_node();
          node->set_op_type("Asinh");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_ACOSH: {
          node = graph->add_node();
          node->set_op_type("Acosh");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_ATANH: {
          node = graph->add_node();
          node->set_op_type("Atanh");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        // ========== Power and Rounding ==========
        case OP_POW: {
          node = graph->add_node();
          node->set_op_type("Pow");
          node->add_input(work_to_onnx[i[0]]);
          node->add_input(work_to_onnx[i[1]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_FABS: {
          node = graph->add_node();
          node->set_op_type("Abs");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_CEIL: {
          node = graph->add_node();
          node->set_op_type("Ceil");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_FLOOR: {
          node = graph->add_node();
          node->set_op_type("Floor");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_SIGN: {
          node = graph->add_node();
          node->set_op_type("Sign");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        // ========== Other Mathematical ==========
        case OP_ERF: {
          node = graph->add_node();
          node->set_op_type("Erf");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        // ========== Tensor Operations ==========
        case OP_MTIMES: {
          node = graph->add_node();
          node->set_op_type("MatMul");
          node->add_input(work_to_onnx[i[0]]);
          node->add_input(work_to_onnx[i[1]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }
        case OP_TRANSPOSE: {
          node = graph->add_node();
          node->set_op_type("Transpose");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_RESHAPE: {
          // Get the MX node to extract target shape info
          MX mx_reshape = f.instruction_MX(k);
          auto sp = mx_reshape.sparsity();

          node = graph->add_node();
          node->set_op_type("Reshape");
          node->add_input(work_to_onnx[i[0]]);

          // Create shape constant as second input
          std::string shape_name = "shape_" + std::to_string(k);
          onnx::NodeProto* shape_node = graph->add_node();
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
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_HORZCAT: {
          node = graph->add_node();
          node->set_op_type("Concat");
          // Add all inputs
          for (casadi_int idx : i) {
            node->add_input(work_to_onnx[idx]);
          }
          node->add_output(node_output);
          // Set axis attribute: axis=1 for horizontal concatenation (columns)
          onnx::AttributeProto* axis_attr = node->add_attribute();
          axis_attr->set_name("axis");
          axis_attr->set_i(1);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_VERTCAT: {
          node = graph->add_node();
          node->set_op_type("Concat");
          // Add all inputs
          for (casadi_int idx : i) {
            node->add_input(work_to_onnx[idx]);
          }
          node->add_output(node_output);
          // Set axis attribute: axis=0 for vertical concatenation (rows)
          onnx::AttributeProto* axis_attr = node->add_attribute();
          axis_attr->set_name("axis");
          axis_attr->set_i(0);
          work_to_onnx[o[0]] = node_output;
          break;
        }

        case OP_HORZSPLIT: {
          // Get the MX node to extract split sizes
          MX mx_split = f.instruction_MX(k);

          node = graph->add_node();
          node->set_op_type("Split");
          node->add_input(work_to_onnx[i[0]]);

          // Set axis attribute: axis=1 for horizontal split (columns)
          onnx::AttributeProto* axis_attr = node->add_attribute();
          axis_attr->set_name("axis");
          axis_attr->set_i(1);

          // Determine split sizes from the output sparsity patterns
          std::vector<casadi_int> split_sizes;
          for (casadi_int j = 0; j < o.size(); ++j) {
            // Get the instruction that produces this output
            MX out_mx = f.instruction_MX(k);
            // For split operations, each output has a size
            auto out_sp = f.sparsity_out(o[j]);
            split_sizes.push_back(out_sp.size2());  // width for horzsplit
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
          for (casadi_int j = 0; j < o.size(); ++j) {
            std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
            node->add_output(output_name);
            work_to_onnx[o[j]] = output_name;
          }
          break;
        }

        case OP_VERTSPLIT: {
          // Get the MX node to extract split sizes
          MX mx_split = f.instruction_MX(k);

          node = graph->add_node();
          node->set_op_type("Split");
          node->add_input(work_to_onnx[i[0]]);

          // Set axis attribute: axis=0 for vertical split (rows)
          onnx::AttributeProto* axis_attr = node->add_attribute();
          axis_attr->set_name("axis");
          axis_attr->set_i(0);

          // Determine split sizes from the output sparsity patterns
          std::vector<casadi_int> split_sizes;
          for (casadi_int j = 0; j < o.size(); ++j) {
            // Get the instruction that produces this output
            MX out_mx = f.instruction_MX(k);
            // For split operations, each output has a size
            auto out_sp = f.sparsity_out(o[j]);
            split_sizes.push_back(out_sp.size1());  // height for vertsplit
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
          for (casadi_int j = 0; j < o.size(); ++j) {
            std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
            node->add_output(output_name);
            work_to_onnx[o[j]] = output_name;
          }
          break;
        }

        case OP_SUBREF: {
          // SUBREF is used for both Slice and GatherElements
          // For now, implement basic slicing support
          // TODO: Distinguish between simple slicing and advanced indexing

          // Basic implementation: assume it's a simple slice operation
          // This is a simplified version - full implementation would need to
          // analyze the slice parameters from the MX node

          if (verbose_) {
            uout() << "ONNX export: OP_SUBREF detected at instruction " << k
                   << " - using simplified Slice export" << std::endl;
          }

          // For now, just pass through as identity
          // A full implementation would create proper Slice nodes
          node = graph->add_node();
          node->set_op_type("Identity");
          node->add_input(work_to_onnx[i[0]]);
          node->add_output(node_output);
          work_to_onnx[o[0]] = node_output;

          casadi_warning("ONNX export: OP_SUBREF (Slice/GatherElements) support is incomplete. "
                        "Using Identity as placeholder.");
          break;
        }

        default:
          casadi_warning("ONNX export: unsupported operation " +
                      std::to_string(op) + " at instruction " + std::to_string(k));
      }
    }

    has_model_ = true;

    if (verbose_) {
      uout() << "Converted CasADi Function to ONNX model: " << f.name() << std::endl;
      uout() << "  Instructions processed: " << n_instr << std::endl;
      uout() << "  ONNX nodes created: " << graph->node_size() << std::endl;
    }
  }

} // namespace casadi
/// \endcond
