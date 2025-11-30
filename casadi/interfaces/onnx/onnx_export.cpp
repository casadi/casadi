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

  // Helper function to add graph inputs
  void add_graph_inputs(onnx::GraphProto* graph, const Function& f) {
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
  }

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
    add_graph_inputs(graph, f);

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

      // Try to process with operation handler
      if (process_operation(graph, f, op, k, i, o, work_to_onnx, node_output)) {
        // Operation was handled successfully
        continue;
      }

      // Handle operations not in common processor
      if (op == OP_CALL) {
          // Function call - could be control flow (If, Loop, Scan) or regular function
          MX mx_call = f.instruction_MX(k);
          Function called_func = mx_call.which_function();

        if (verbose_) {
          uout() << "ONNX export: Found function call to '" << called_func.name()
                 << "' at instruction " << k << std::endl;
        }

        // Check if it's a control flow operation
        if (is_if_else_function(called_func)) {
          // TODO: Implement If export
          casadi_error("ONNX export: If/conditional operators not yet fully implemented for export. "
                      "Import is supported but export of control flow requires additional work.");

        } else if (is_mapaccum_function(called_func)) {
          // TODO: Implement Loop export
          casadi_error("ONNX export: Loop/mapaccum operators not yet fully implemented for export. "
                      "Import is supported but export of control flow requires additional work.");

        } else if (is_map_function(called_func)) {
          // TODO: Implement Scan export
          casadi_error("ONNX export: Scan/map operators not yet fully implemented for export. "
                      "Import is supported but export of control flow requires additional work.");

        } else {
          // Regular function call
          casadi_error("ONNX export: Function calls to '" + called_func.name() +
                      "' cannot be exported. ONNX does not support arbitrary function calls. "
                      "Consider inlining the function before export.");
        }
      } else {
        // Unknown/unsupported operation
        casadi_error("ONNX export: unsupported operation code " +
                    std::to_string(op) + " at instruction " + std::to_string(k) +
                    ". The CasADi Function contains operations that cannot be exported to ONNX.");
      }
    }

    has_model_ = true;

    if (verbose_) {
      uout() << "Converted CasADi Function to ONNX model: " << f.name() << std::endl;
      uout() << "  Instructions processed: " << n_instr << std::endl;
      uout() << "  ONNX nodes created: " << graph->node_size() << std::endl;
    }
  }

  // ========== Control Flow Support ==========

  bool OnnxTranslator::is_if_else_function(const Function& f) const {
    // Check if function name contains "if_else" or "conditional"
    std::string fname = f.name();
    return fname.find("if_else") != std::string::npos ||
           fname.find("conditional") != std::string::npos ||
           fname.find("if_") == 0;
  }

  bool OnnxTranslator::is_mapaccum_function(const Function& f) const {
    // Check if function name contains "mapaccum" or "accum"
    std::string fname = f.name();
    return fname.find("mapaccum") != std::string::npos ||
           fname.find("accum") != std::string::npos;
  }

  bool OnnxTranslator::is_map_function(const Function& f) const {
    // Check if function name contains "map" but not "mapaccum"
    std::string fname = f.name();
    return fname.find("map") != std::string::npos &&
           fname.find("mapaccum") == std::string::npos &&
           fname.find("accum") == std::string::npos;
  }

  onnx::GraphProto* OnnxTranslator::function_to_graph(
      const Function& f,
      const std::string& graph_name) {

    // Create a new graph
    onnx::GraphProto* graph = new onnx::GraphProto();
    graph->set_name(graph_name);

    if (verbose_) {
      uout() << "  Converting Function '" << f.name() << "' to ONNX subgraph '"
             << graph_name << "'" << std::endl;
      uout() << "    Inputs: " << f.n_in() << ", Outputs: " << f.n_out()
             << ", Instructions: " << f.n_instructions() << std::endl;
    }

    // Add inputs to subgraph
    add_graph_inputs(graph, f);

    // Map from work vector index to ONNX node name
    std::map<casadi_int, std::string> work_to_onnx;

    // Process instructions
    casadi_int n_instr = f.n_instructions();
    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i_vec = f.instruction_input(k);

      std::string node_output = "n" + std::to_string(k);

      // Try to process with operation handler
      if (process_operation(graph, f, op, k, i_vec, o, work_to_onnx, node_output)) {
        // Operation was handled successfully
        continue;
      }

      // Handle operations not in common processor
      if (op == OP_CALL) {
        // Get the called function
        MX mx_call = f.instruction_MX(k);
        Function called_func = mx_call.which_function();

        if (verbose_) {
          uout() << "    Found OP_CALL: " << called_func.name() << std::endl;
        }

        // For now, error on nested function calls in subgraphs
        // Full implementation would recursively handle this
        casadi_error("ONNX export: Nested function calls in subgraphs not yet supported. "
                    "Found call to '" + called_func.name() + "' in subgraph '" +
                    graph_name + "'.");
      } else {
        // For unsupported operations in subgraphs, provide helpful error
        casadi_error("ONNX export: Unsupported operation code " + std::to_string(op) +
                    " in subgraph '" + graph_name + "' at instruction " +
                    std::to_string(k) + ". Subgraphs have limited operation support.");
      }
    }

    if (verbose_) {
      uout() << "    Created " << graph->node_size() << " nodes in subgraph" << std::endl;
    }

    return graph;
  }

} // namespace casadi
/// \endcond
