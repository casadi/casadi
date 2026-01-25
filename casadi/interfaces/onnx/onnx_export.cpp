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
  void add_graph_inputs(onnx::GraphProto* graph, const Function& f,
                        const std::string& name_prefix) {
    for (casadi_int i = 0; i < f.n_in(); ++i) {
      onnx::ValueInfoProto* input = graph->add_input();
      std::string input_name;

      if (!name_prefix.empty()) {
        // For subgraphs, ALWAYS use prefix to ensure uniqueness
        input_name = name_prefix + "_i" + std::to_string(i);
      } else {
        // For main graph, use function's input name if available
        input_name = f.name_in(i);
        if (input_name.empty()) {
          input_name = "input_" + std::to_string(i);
        }
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

  // Helper function to add graph outputs
  void add_graph_outputs(onnx::GraphProto* graph, const Function& f,
                         const std::string& name_prefix) {
    for (casadi_int i = 0; i < f.n_out(); ++i) {
      onnx::ValueInfoProto* output = graph->add_output();
      std::string output_name;

      if (!name_prefix.empty()) {
        // For subgraphs, ALWAYS use prefix to ensure uniqueness
        // (don't use function's original output names as they might conflict)
        output_name = name_prefix + "_o" + std::to_string(i);
      } else {
        // For main graph, use function's output name if available
        output_name = f.name_out(i);
        if (output_name.empty()) {
          output_name = "output_" + std::to_string(i);
        }
      }
      output->set_name(output_name);

      // Set tensor type and shape
      onnx::TypeProto* type = output->mutable_type();
      onnx::TypeProto::Tensor* tensor_type = type->mutable_tensor_type();
      tensor_type->set_elem_type(onnx::TensorProto::DOUBLE);

      // Add shape dimensions
      onnx::TensorShapeProto* shape = tensor_type->mutable_shape();
      auto sp = f.sparsity_out(i);
      shape->add_dim()->set_dim_value(sp.size1());
      shape->add_dim()->set_dim_value(sp.size2());
    }
  }

  void OnnxTranslator::load(const Function& f) {
    // Create ONNX model from CasADi Function
    model_.Clear();
    exported_functions_.clear();  // Clear previously exported functions
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

    // Pre-scan: identify outputs with multiple segments (for vertcat outputs)
    // Map: output_idx -> count of segments
    std::map<casadi_int, casadi_int> output_segment_count;
    for (casadi_int k = 0; k < n_instr; ++k) {
      if (f.instruction_id(k) == OP_OUTPUT) {
        MX mx = f.instruction_MX(k);
        Dict info = mx.info();
        casadi_int output_idx = info["ind"];
        output_segment_count[output_idx]++;
      }
    }

    // Track segment values for multi-segment outputs: output_idx -> (offset -> onnx_name)
    std::map<casadi_int, std::map<casadi_int, std::string>> output_segment_values;

    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i = f.instruction_input(k);

      // Node name for this operation's result - use instruction index for uniqueness
      std::string node_output = "n" + std::to_string(k);

      // Handle OP_OUTPUT specially to support multi-segment outputs (vertcat)
      if (op == OP_OUTPUT) {
        MX mx = f.instruction_MX(k);
        Dict info = mx.info();
        casadi_int output_idx = info["ind"];
        casadi_int offset = info["offset"];

        std::string output_name = f.name_out(output_idx);
        if (output_name.empty()) {
          output_name = "output_" + std::to_string(output_idx);
        }

        std::string input_onnx_name = work_to_onnx[i[0]];

        if (output_segment_count[output_idx] > 1) {
          // Multi-segment output: write to temp name, will Concat later
          std::string seg_name = output_name + "_seg" + std::to_string(offset);
          onnx::NodeProto* id_node = graph->add_node();
          id_node->set_op_type("Identity");
          id_node->add_input(input_onnx_name);
          id_node->add_output(seg_name);
          output_segment_values[output_idx][offset] = seg_name;
        } else {
          // Single-segment output: write directly to output name
          onnx::NodeProto* id_node = graph->add_node();
          id_node->set_op_type("Identity");
          id_node->add_input(input_onnx_name);
          id_node->add_output(output_name);
        }
        continue;
      }

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

        // Check if it's a control flow operation (not supported)
        if (is_if_else_function(called_func)) {
          casadi_error("ONNX export: Function.if_else is not supported.");
        } else if (is_mapaccum_function(called_func)) {
          casadi_error("ONNX export: mapaccum (Loop) is not supported.");
        } else if (is_map_function(called_func)) {
          casadi_error("ONNX export: Function.map (Scan) is not supported.");
        } else {
          // Regular function call - export as ONNX local function using FunctionProto
          std::string func_name = called_func.name();
          std::string domain = "casadi";

          if (verbose_) {
            uout() << "  Exporting regular function call to '" << func_name
                   << "' (domain: " << domain << ")" << std::endl;
          }

          // Export function if not already done
          if (exported_functions_.find(func_name) == exported_functions_.end()) {
            onnx::FunctionProto* func_proto = function_to_function_proto(called_func, domain);
            *model_.add_functions() = *func_proto;
            delete func_proto;
            exported_functions_.insert(func_name);
          }

          // Create node that calls the function
          onnx::NodeProto* call_node = graph->add_node();
          call_node->set_op_type(func_name);
          call_node->set_domain(domain);

          // Connect inputs
          for (size_t j = 0; j < i.size(); ++j) {
            call_node->add_input(work_to_onnx[i[j]]);
          }

          // Connect outputs
          for (size_t j = 0; j < o.size(); ++j) {
            std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
            call_node->add_output(output_name);
            work_to_onnx[o[j]] = output_name;
          }

          if (verbose_) {
            uout() << "    Created function call node with " << o.size() << " outputs" << std::endl;
          }

          continue;
        }
      } else {
        // Unknown/unsupported operation
        casadi_error("ONNX export: unsupported operation code " +
                    std::to_string(op) + " at instruction " + std::to_string(k) +
                    ". The CasADi Function contains operations that cannot be exported to ONNX.");
      }
    }

    // Create Concat nodes for multi-segment outputs (vertcat outputs)
    for (const auto& kv : output_segment_values) {
      casadi_int output_idx = kv.first;
      const auto& segments = kv.second;  // map<offset, onnx_name>, already sorted by offset

      std::string output_name = f.name_out(output_idx);
      if (output_name.empty()) {
        output_name = "output_" + std::to_string(output_idx);
      }

      // Create Concat node to combine segments
      onnx::NodeProto* concat_node = graph->add_node();
      concat_node->set_op_type("Concat");

      // Add inputs in order of offset (map is sorted by key)
      for (const auto& seg : segments) {
        concat_node->add_input(seg.second);
      }

      concat_node->add_output(output_name);

      // Set axis=0 for vertical concatenation
      onnx::AttributeProto* axis_attr = concat_node->add_attribute();
      axis_attr->set_name("axis");
      axis_attr->set_type(onnx::AttributeProto::INT);
      axis_attr->set_i(0);
    }

    // Add graph outputs (for main graph, use empty prefix to use function's output names)
    add_graph_outputs(graph, f, "");

    // If we exported any functions, add the casadi domain opset_import
    if (!exported_functions_.empty()) {
      onnx::OperatorSetIdProto* casadi_opset = model_.add_opset_import();
      casadi_opset->set_domain("casadi");
      casadi_opset->set_version(1);

      if (verbose_) {
        uout() << "  Exported " << exported_functions_.size()
               << " function(s) to casadi domain" << std::endl;
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
    // Check using class_name for reliable detection
    return f.class_name() == "Switch";
  }

  bool OnnxTranslator::is_mapaccum_function(const Function& f) const {
    // Both mapaccum and map have class_name "Map" or "OmpMap"
    // Use name heuristic to differentiate
    std::string class_name = f.class_name();
    if (class_name != "Map" && class_name != "OmpMap") {
      return false;
    }

    // Check function name for "mapaccum" or "accum"
    std::string fname = f.name();
    return fname.find("mapaccum") != std::string::npos ||
           fname.find("accum") != std::string::npos;
  }

  bool OnnxTranslator::is_map_function(const Function& f) const {
    // Both mapaccum and map have class_name "Map" or "OmpMap"
    // This is map if it's Map class but NOT mapaccum
    std::string class_name = f.class_name();
    if (class_name != "Map" && class_name != "OmpMap") {
      return false;
    }

    // Exclude mapaccum by checking name doesn't contain "mapaccum" or "accum"
    std::string fname = f.name();
    return fname.find("mapaccum") == std::string::npos &&
           fname.find("accum") == std::string::npos;
  }

  onnx::FunctionProto* OnnxTranslator::function_to_function_proto(
      const Function& f,
      const std::string& domain) {

    // Create a new FunctionProto
    onnx::FunctionProto* func = new onnx::FunctionProto();
    func->set_name(f.name());
    func->set_domain(domain);

    if (verbose_) {
      uout() << "  Converting Function '" << f.name() << "' to ONNX FunctionProto"
             << " (domain: " << domain << ")" << std::endl;
      uout() << "    Inputs: " << f.n_in() << ", Outputs: " << f.n_out()
             << ", Instructions: " << f.n_instructions() << std::endl;
    }

    // Add function inputs (parameter names)
    for (casadi_int i = 0; i < f.n_in(); ++i) {
      std::string input_name = f.name_in(i);
      if (input_name.empty()) {
        input_name = "input_" + std::to_string(i);
      }
      func->add_input(input_name);
    }

    // Add function outputs (parameter names)
    for (casadi_int i = 0; i < f.n_out(); ++i) {
      std::string output_name = f.name_out(i);
      if (output_name.empty()) {
        output_name = "output_" + std::to_string(i);
      }
      func->add_output(output_name);
    }

    // Map from work vector index to ONNX node name
    std::map<casadi_int, std::string> work_to_onnx;

    // Process instructions and add nodes to the function
    casadi_int n_instr = f.n_instructions();
    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i_vec = f.instruction_input(k);

      std::string node_output = "n" + std::to_string(k);

      // Handle OP_INPUT specially for FunctionProto - just map to parameter name (no Identity node)
      if (op == OP_INPUT) {
        casadi_int input_idx = i_vec[0];
        std::string input_name = f.name_in(input_idx);
        if (input_name.empty()) {
          input_name = "input_" + std::to_string(input_idx);
        }
        work_to_onnx[o[0]] = input_name;
        continue;
      }

      // Handle OP_OUTPUT specially for FunctionProto - use declared output names
      if (op == OP_OUTPUT) {
        casadi_int output_idx = o[0];
        std::string output_name = f.name_out(output_idx);
        if (output_name.empty()) {
          output_name = "output_" + std::to_string(output_idx);
        }

        // Create Identity node to connect internal result to output name
        onnx::NodeProto* id_node = func->add_node();
        id_node->set_op_type("Identity");
        id_node->add_input(work_to_onnx[i_vec[0]]);
        id_node->add_output(output_name);
        continue;
      }

      // Handle OP_CALL - recursive function call
      if (op == OP_CALL) {
        MX mx_call = f.instruction_MX(k);
        Function called_func = mx_call.which_function();

        if (verbose_) {
          uout() << "    Found OP_CALL: " << called_func.name() << std::endl;
        }

        // Check for control flow - these need special handling
        if (is_if_else_function(called_func) ||
            is_mapaccum_function(called_func) ||
            is_map_function(called_func)) {
          casadi_error("ONNX export: Control flow in nested functions not yet supported. "
                      "Found '" + called_func.name() + "' in function '" + f.name() + "'.");
        }

        // Regular function call - export the called function if not already done
        std::string func_name = called_func.name();
        if (exported_functions_.find(func_name) == exported_functions_.end()) {
          // Recursively export the called function
          onnx::FunctionProto* nested_func = function_to_function_proto(called_func, domain);
          *model_.add_functions() = *nested_func;
          delete nested_func;
          exported_functions_.insert(func_name);
        }

        // Create node that calls the function
        onnx::NodeProto* call_node = func->add_node();
        call_node->set_op_type(func_name);
        call_node->set_domain(domain);

        // Connect inputs
        for (size_t j = 0; j < i_vec.size(); ++j) {
          call_node->add_input(work_to_onnx[i_vec[j]]);
        }

        // Connect outputs
        for (size_t j = 0; j < o.size(); ++j) {
          std::string output = "call_" + func_name + "_" + std::to_string(k) + "_out" + std::to_string(j);
          call_node->add_output(output);
          work_to_onnx[o[j]] = output;
        }

        continue;
      }

      // All other operations use shared implementation via callback
      if (process_operation([&]() { return func->add_node(); },
                            f, op, k, i_vec, o, work_to_onnx, node_output)) {
        continue;
      }

      // Unsupported operation
      casadi_error("ONNX export: Unsupported operation code " + std::to_string(op) +
                  " in function '" + f.name() + "' at instruction " + std::to_string(k));
    }

    if (verbose_) {
      uout() << "    Created " << func->node_size() << " nodes in FunctionProto" << std::endl;
    }

    return func;
  }

} // namespace casadi
/// \endcond
