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
          // If/else export
          if (verbose_) {
            uout() << "  Exporting If/else operator" << std::endl;
          }

          // Extract branches from Switch using info()
          Dict info = called_func.info();
          Function f_true = info["f_def"];                    // Default/true branch
          std::vector<Function> f_cases = info["f"];          // Cases vector
          Function f_false = f_cases[0];                      // First case is false branch

          // For If branches, we pass the branch input values as outer scope variables
          // The Switch function has inputs: [condition, value1, value2, ...]
          // The branches need access to the values (not the condition)
          std::vector<std::string> outer_scope_vars;
          for (size_t j = 1; j < i.size(); ++j) {
            outer_scope_vars.push_back(work_to_onnx[i[j]]);
          }

          if (verbose_) {
            uout() << "    Branch outer scope vars: ";
            for (const auto& v : outer_scope_vars) {
              uout() << v << " ";
            }
            uout() << std::endl;
          }

          // Convert branches to ONNX subgraphs with outer scope variable references
          onnx::GraphProto* then_graph = function_to_graph(f_true, "then_branch", outer_scope_vars);
          onnx::GraphProto* else_graph = function_to_graph(f_false, "else_branch", outer_scope_vars);

          // Create ONNX If node
          onnx::NodeProto* if_node = graph->add_node();
          if_node->set_op_type("If");

          // Add condition as input (first input to the Switch function call)
          if_node->add_input(work_to_onnx[i[0]]);

          // Add branches as graph attributes
          onnx::AttributeProto* then_attr = if_node->add_attribute();
          then_attr->set_name("then_branch");
          then_attr->set_type(onnx::AttributeProto::GRAPH);
          then_attr->set_allocated_g(then_graph);

          onnx::AttributeProto* else_attr = if_node->add_attribute();
          else_attr->set_name("else_branch");
          else_attr->set_type(onnx::AttributeProto::GRAPH);
          else_attr->set_allocated_g(else_graph);

          // Add outputs from the If node
          for (casadi_int j = 0; j < o.size(); ++j) {
            std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
            if_node->add_output(output_name);
            work_to_onnx[o[j]] = output_name;
          }

          if (verbose_) {
            uout() << "  Created If node with " << o.size() << " outputs" << std::endl;
          }

          continue;

        } else if (is_mapaccum_function(called_func)) {
          // Loop/mapaccum export
          if (verbose_) {
            uout() << "  Exporting Loop/mapaccum operator" << std::endl;
          }

          // Extract base function from Map using get_function
          Function base_func = called_func.get_function("f");
          Dict info = called_func.info();
          casadi_int n_iter = info["n"];  // Number of iterations

          if (verbose_) {
            uout() << "    Base function: " << base_func.name()
                   << ", iterations: " << n_iter << std::endl;
          }

          // Create max_iter constant
          std::string max_iter_name = "max_iter_" + std::to_string(k);
          onnx::NodeProto* max_iter_node = graph->add_node();
          max_iter_node->set_op_type("Constant");
          max_iter_node->add_output(max_iter_name);

          onnx::AttributeProto* max_iter_attr = max_iter_node->add_attribute();
          max_iter_attr->set_name("value");
          onnx::TensorProto* max_iter_tensor = max_iter_attr->mutable_t();
          max_iter_tensor->set_data_type(onnx::TensorProto::INT64);
          max_iter_tensor->add_dims(1);  // Scalar
          max_iter_tensor->add_int64_data(n_iter);

          // Create initial condition (always true for fixed iteration count)
          std::string init_cond_name = "init_cond_" + std::to_string(k);
          onnx::NodeProto* init_cond_node = graph->add_node();
          init_cond_node->set_op_type("Constant");
          init_cond_node->add_output(init_cond_name);

          onnx::AttributeProto* init_cond_attr = init_cond_node->add_attribute();
          init_cond_attr->set_name("value");
          onnx::TensorProto* init_cond_tensor = init_cond_attr->mutable_t();
          init_cond_tensor->set_data_type(onnx::TensorProto::BOOL);
          init_cond_tensor->add_dims(1);  // Scalar
          init_cond_tensor->add_int32_data(1);  // true

          // Convert body to ONNX subgraph
          // NOTE: CasADi mapaccum body may not match ONNX Loop signature exactly
          // ONNX expects: [iter_num, cond, state...] -> [cond_out, state_out...]
          // For now, export as-is and rely on compatible structure
          onnx::GraphProto* body_graph = function_to_graph(base_func, "loop_body");

          // Create ONNX Loop node
          onnx::NodeProto* loop_node = graph->add_node();
          loop_node->set_op_type("Loop");

          // Add max_iter as first input
          loop_node->add_input(max_iter_name);

          // Add initial condition as second input
          loop_node->add_input(init_cond_name);

          // Add loop-carried dependencies as remaining inputs
          for (casadi_int j = 0; j < i.size(); ++j) {
            loop_node->add_input(work_to_onnx[i[j]]);
          }

          // Add body as graph attribute
          onnx::AttributeProto* body_attr = loop_node->add_attribute();
          body_attr->set_name("body");
          body_attr->set_type(onnx::AttributeProto::GRAPH);
          body_attr->set_allocated_g(body_graph);

          // Add outputs from the Loop node
          for (casadi_int j = 0; j < o.size(); ++j) {
            std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
            loop_node->add_output(output_name);
            work_to_onnx[o[j]] = output_name;
          }

          if (verbose_) {
            uout() << "  Created Loop node with " << o.size() << " outputs" << std::endl;
          }

          continue;

        } else if (is_map_function(called_func)) {
          // Scan/map export
          if (verbose_) {
            uout() << "  Exporting Scan/map operator" << std::endl;
          }

          // Extract base function from Map using get_function
          Function base_func = called_func.get_function("f");
          Dict info = called_func.info();
          casadi_int n_iter = info["n"];  // Number of iterations

          // Determine num_scan_inputs by analyzing input sparsity
          // Scan inputs have dimension matching n_iter
          casadi_int n_base_inputs = base_func.n_in();
          casadi_int num_scan_inputs = 0;

          // Count inputs from the end that are scan inputs (size matches n_iter)
          for (casadi_int j = n_base_inputs - 1; j >= 0; --j) {
            auto sp = base_func.sparsity_in(j);
            // Scan inputs are horizontally concatenated (cols == n_iter)
            if (sp.size2() == n_iter) {
              num_scan_inputs++;
            } else {
              break;  // Found first state input, stop counting
            }
          }

          casadi_int num_state_inputs = n_base_inputs - num_scan_inputs;

          if (verbose_) {
            uout() << "    Base function: " << base_func.name()
                   << ", state inputs: " << num_state_inputs
                   << ", scan inputs: " << num_scan_inputs << std::endl;
          }

          // Convert body to ONNX subgraph
          onnx::GraphProto* body_graph = function_to_graph(base_func, "scan_body");

          // Create ONNX Scan node
          onnx::NodeProto* scan_node = graph->add_node();
          scan_node->set_op_type("Scan");

          // Add inputs (state variables + scan inputs)
          for (casadi_int j = 0; j < i.size(); ++j) {
            scan_node->add_input(work_to_onnx[i[j]]);
          }

          // Add body as graph attribute
          onnx::AttributeProto* body_attr = scan_node->add_attribute();
          body_attr->set_name("body");
          body_attr->set_type(onnx::AttributeProto::GRAPH);
          body_attr->set_allocated_g(body_graph);

          // Add num_scan_inputs attribute
          onnx::AttributeProto* num_scan_attr = scan_node->add_attribute();
          num_scan_attr->set_name("num_scan_inputs");
          num_scan_attr->set_type(onnx::AttributeProto::INT);
          num_scan_attr->set_i(num_scan_inputs);

          // Add outputs from the Scan node
          for (casadi_int j = 0; j < o.size(); ++j) {
            std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
            scan_node->add_output(output_name);
            work_to_onnx[o[j]] = output_name;
          }

          if (verbose_) {
            uout() << "  Created Scan node with " << o.size() << " outputs" << std::endl;
          }

          continue;

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

  onnx::GraphProto* OnnxTranslator::function_to_graph(
      const Function& f,
      const std::string& graph_name,
      const std::vector<std::string>& outer_scope_inputs) {

    // Create a new graph
    onnx::GraphProto* graph = new onnx::GraphProto();
    graph->set_name(graph_name);

    if (verbose_) {
      uout() << "  Converting Function '" << f.name() << "' to ONNX subgraph '"
             << graph_name << "'" << std::endl;
      uout() << "    Inputs: " << f.n_in() << ", Outputs: " << f.n_out()
             << ", Instructions: " << f.n_instructions() << std::endl;
      if (!outer_scope_inputs.empty()) {
        uout() << "    Outer scope inputs: ";
        for (const auto& name : outer_scope_inputs) {
          uout() << name << " ";
        }
        uout() << std::endl;
      }
    }

    // For subgraphs with outer scope inputs (like If branches), don't add graph inputs
    // The subgraph will reference outer scope variables directly
    if (outer_scope_inputs.empty()) {
      add_graph_inputs(graph, f, "");
    }

    // Map from work vector index to ONNX node name
    std::map<casadi_int, std::string> work_to_onnx;

    // If we have outer scope inputs, pre-populate work_to_onnx so OP_INPUT
    // instructions will find the correct outer scope variable names
    for (size_t i = 0; i < outer_scope_inputs.size(); ++i) {
      // Map the function's input work index to the outer scope variable name
      // For a function with n inputs, inputs are at work indices 0 to n-1
      work_to_onnx[static_cast<casadi_int>(i)] = outer_scope_inputs[i];
    }

    // Process instructions
    casadi_int n_instr = f.n_instructions();
    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i_vec = f.instruction_input(k);

      std::string node_output = "n" + std::to_string(k);

      // Special handling for OP_INPUT when using outer scope inputs
      // (for If branch subgraphs that reference outer scope variables)
      if (op == OP_INPUT && !outer_scope_inputs.empty()) {
        casadi_int input_idx = i_vec[0];  // Function input index
        if (static_cast<size_t>(input_idx) < outer_scope_inputs.size()) {
          // Map the output work index directly to the outer scope variable name
          // No need to create an Identity node - we just reference the outer variable
          work_to_onnx[o[0]] = outer_scope_inputs[input_idx];
          continue;
        }
      }

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
        std::string func_name = called_func.name();
        std::string domain = "casadi";

        if (verbose_) {
          uout() << "    Found OP_CALL: " << func_name << " in subgraph" << std::endl;
        }

        // Check for control flow - these need special handling within subgraphs
        if (is_if_else_function(called_func) ||
            is_mapaccum_function(called_func) ||
            is_map_function(called_func)) {
          casadi_error("ONNX export: Control flow in subgraphs not yet supported. "
                      "Found '" + func_name + "' in subgraph '" + graph_name + "'.");
        }

        // Regular function call - export as ONNX local function
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
        for (size_t j = 0; j < i_vec.size(); ++j) {
          call_node->add_input(work_to_onnx[i_vec[j]]);
        }

        // Connect outputs
        for (size_t j = 0; j < o.size(); ++j) {
          std::string output = node_output + "_out" + std::to_string(j);
          call_node->add_output(output);
          work_to_onnx[o[j]] = output;
        }

        continue;
      } else {
        // For unsupported operations in subgraphs, provide helpful error
        casadi_error("ONNX export: Unsupported operation code " + std::to_string(op) +
                    " in subgraph '" + graph_name + "' at instruction " +
                    std::to_string(k) + ". Subgraphs have limited operation support.");
      }
    }

    // Add graph outputs (without prefix - subgraphs are self-contained)
    add_graph_outputs(graph, f, "");

    if (verbose_) {
      uout() << "    Created " << graph->node_size() << " nodes in subgraph" << std::endl;
    }

    return graph;
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

      // Handle OP_CONST - create constant node
      if (op == OP_CONST) {
        MX mx_const = f.instruction_MX(k);
        DM dm_const = static_cast<DM>(mx_const);

        onnx::NodeProto* const_node = func->add_node();
        const_node->set_op_type("Constant");
        const_node->add_output(node_output);

        onnx::AttributeProto* attr = const_node->add_attribute();
        attr->set_name("value");
        onnx::TensorProto* tensor = attr->mutable_t();
        tensor->set_data_type(onnx::TensorProto::DOUBLE);

        tensor->add_dims(dm_const.size1());
        if (dm_const.size2() > 1) tensor->add_dims(dm_const.size2());

        for (casadi_int idx = 0; idx < dm_const.numel(); ++idx) {
          tensor->add_double_data(static_cast<double>(dm_const->at(idx)));
        }

        work_to_onnx[o[0]] = node_output;
        continue;
      }

      // Try simple operations using centralized lookup table
      const OpMapping* mapping = get_op_mapping(op);
      if (mapping) {
        onnx::NodeProto* node = func->add_node();
        node->set_op_type(mapping->onnx_name);
        if (mapping->arity == 1 && i_vec.size() >= 1) {
          node->add_input(work_to_onnx[i_vec[0]]);
        } else if (mapping->arity == 2 && i_vec.size() >= 2) {
          node->add_input(work_to_onnx[i_vec[0]]);
          node->add_input(work_to_onnx[i_vec[1]]);
        }
        node->add_output(node_output);
        work_to_onnx[o[0]] = node_output;
        continue;
      }

      // Handle OP_SQ (square): x^2 = x * x
      if (op == OP_SQ && i_vec.size() == 1 && o.size() == 1) {
        onnx::NodeProto* node = func->add_node();
        node->set_op_type("Mul");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_output(node_output);
        work_to_onnx[o[0]] = node_output;
        continue;
      }

      // Handle OP_TWICE (double): 2*x
      if (op == OP_TWICE && i_vec.size() == 1 && o.size() == 1) {
        std::string const_name = "const_2_" + std::to_string(k);
        onnx::NodeProto* const_node = func->add_node();
        const_node->set_op_type("Constant");
        const_node->add_output(const_name);
        onnx::AttributeProto* attr = const_node->add_attribute();
        attr->set_name("value");
        onnx::TensorProto* tensor = attr->mutable_t();
        tensor->set_data_type(onnx::TensorProto::DOUBLE);
        tensor->add_double_data(2.0);

        onnx::NodeProto* mul_node = func->add_node();
        mul_node->set_op_type("Mul");
        mul_node->add_input(const_name);
        mul_node->add_input(work_to_onnx[i_vec[0]]);
        mul_node->add_output(node_output);
        work_to_onnx[o[0]] = node_output;
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
