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

  std::vector<std::string> OnnxTranslator::analyze_outer_scope_dependencies(
      const onnx::GraphProto& subgraph,
      const std::set<std::string>& available_in_scope) const {

    // Step 1: Collect all locally-defined names in the subgraph
    std::set<std::string> local_names;

    // Subgraph inputs are local
    for (int i = 0; i < subgraph.input_size(); ++i) {
      local_names.insert(subgraph.input(i).name());
    }

    // Subgraph initializers are local
    for (int i = 0; i < subgraph.initializer_size(); ++i) {
      local_names.insert(subgraph.initializer(i).name());
    }

    // Node outputs are local
    for (int i = 0; i < subgraph.node_size(); ++i) {
      const onnx::NodeProto& node = subgraph.node(i);
      for (int j = 0; j < node.output_size(); ++j) {
        if (!node.output(j).empty()) {
          local_names.insert(node.output(j));
        }
      }
    }

    // Step 2: Find outer scope dependencies
    std::set<std::string> outer_deps;

    for (int i = 0; i < subgraph.node_size(); ++i) {
      const onnx::NodeProto& node = subgraph.node(i);

      // Check each input of this node
      for (int j = 0; j < node.input_size(); ++j) {
        const std::string& input_name = node.input(j);

        // Skip empty inputs (optional parameters in ONNX)
        if (input_name.empty()) continue;

        // If not defined locally, must come from outer scope
        if (local_names.count(input_name) == 0) {
          // Validate it exists in outer scope
          if (available_in_scope.count(input_name)) {
            outer_deps.insert(input_name);
          } else {
            casadi_error("Subgraph references undefined variable '" + input_name +
                       "' in node '" + node.name() + "' (op: " + node.op_type() + ")");
          }
        }
      }

      // Step 3: Recursively handle nested subgraphs
      for (int a = 0; a < node.attribute_size(); ++a) {
        const onnx::AttributeProto& attr = node.attribute(a);

        if (attr.type() == onnx::AttributeProto::GRAPH) {
          // Nested subgraph can see current subgraph's scope
          std::set<std::string> extended_scope = available_in_scope;
          extended_scope.insert(local_names.begin(), local_names.end());

          // Recursively analyze nested subgraph
          std::vector<std::string> nested_deps =
              analyze_outer_scope_dependencies(attr.g(), extended_scope);

          // Nested deps that aren't local to us become our outer deps
          for (const auto& dep : nested_deps) {
            if (local_names.count(dep) == 0) {
              outer_deps.insert(dep);
            }
          }
        }
      }
    }

    // Convert set to vector and return
    return std::vector<std::string>(outer_deps.begin(), outer_deps.end());
  }

  Function OnnxTranslator::translate_subgraph_to_function(
      const onnx::GraphProto& subgraph,
      const std::string& function_name,
      const std::map<std::string, MX>& outer_scope_vars,
      const std::vector<std::string>& outer_deps) {

    // Create local value_map for subgraph scope
    std::map<std::string, MX> value_map;

    // Step 1: Add outer scope dependencies to value_map
    for (const auto& dep_name : outer_deps) {
      auto it = outer_scope_vars.find(dep_name);
      casadi_assert(it != outer_scope_vars.end(),
                   "Outer scope dependency '" + dep_name + "' not found");
      value_map[dep_name] = it->second;
    }

    // Step 2: Process subgraph initializers (constants)
    for (int i = 0; i < subgraph.initializer_size(); ++i) {
      const onnx::TensorProto& tensor = subgraph.initializer(i);
      std::string tensor_name = tensor.name();
      DM dm_const = tensor_to_dm(tensor);
      value_map[tensor_name] = MX(dm_const);
    }

    // Step 3: Create MX symbols for subgraph inputs
    std::vector<MX> func_inputs;
    std::vector<std::string> input_names;

    for (int i = 0; i < subgraph.input_size(); ++i) {
      const onnx::ValueInfoProto& input = subgraph.input(i);
      std::string input_name = input.name();

      // Skip if already in value_map (it's an initializer or outer scope var)
      if (value_map.count(input_name)) {
        continue;
      }

      // Extract shape
      const onnx::TensorShapeProto& shape = input.type().tensor_type().shape();
      casadi_int rows = get_dimension(shape, 0);
      casadi_int cols = get_dimension(shape, 1);

      // Create MX symbol
      MX mx_input = MX::sym(input_name, rows, cols);
      value_map[input_name] = mx_input;
      func_inputs.push_back(mx_input);
      input_names.push_back(input_name);
    }

    // Add outer scope dependencies as function inputs
    for (const auto& dep_name : outer_deps) {
      func_inputs.push_back(value_map[dep_name]);
      input_names.push_back(dep_name);
    }

    // Step 4: Process subgraph nodes (using same logic as create())
    for (int i = 0; i < subgraph.node_size(); ++i) {
      const onnx::NodeProto& node = subgraph.node(i);
      std::string op_type = node.op_type();

      // Gather input tensors
      std::vector<MX> node_inputs;
      for (int j = 0; j < node.input_size(); ++j) {
        std::string input_name = node.input(j);

        // Empty string means optional input not provided
        if (input_name.empty()) {
          node_inputs.push_back(MX());
          continue;
        }

        casadi_assert(value_map.count(input_name),
                     "Unknown input tensor '" + input_name +
                     "' in subgraph '" + function_name + "' node " +
                     std::to_string(i) + " (op_type: " + op_type + ")");
        node_inputs.push_back(value_map[input_name]);
      }

      // Compute output using shared operation logic
      MX output = process_node_operation(op_type, node, node_inputs);

      // Store output
      casadi_assert(node.output_size() >= 1,
                   "Node must have at least one output");
      std::string output_name = node.output(0);
      value_map[output_name] = output;
    }

    // Step 5: Extract subgraph outputs
    std::vector<MX> func_outputs;
    std::vector<std::string> output_names;

    for (int i = 0; i < subgraph.output_size(); ++i) {
      const onnx::ValueInfoProto& output_info = subgraph.output(i);
      std::string output_name = output_info.name();

      casadi_assert(value_map.count(output_name),
                   "Subgraph output '" + output_name + "' not found in value_map");

      func_outputs.push_back(value_map[output_name]);
      output_names.push_back(output_name);
    }

    // Step 6: Create and return Function
    return Function(function_name, func_inputs, func_outputs,
                   input_names, output_names);
  }

  Function OnnxTranslator::translate_loop_body_to_function(
      const onnx::GraphProto& body,
      const std::string& function_name,
      const std::map<std::string, MX>& outer_scope_vars,
      const std::vector<std::string>& outer_deps) {

    // Loop body has special structure:
    // Inputs: [iteration_num, cond_in, loop_carried_deps...]
    // Outputs: [cond_out, loop_carried_deps..., scan_outputs...]

    // Create local value_map for body scope
    std::map<std::string, MX> value_map;

    // Step 1: Add outer scope dependencies to value_map
    for (const auto& dep_name : outer_deps) {
      auto it = outer_scope_vars.find(dep_name);
      casadi_assert(it != outer_scope_vars.end(),
                   "Outer scope dependency '" + dep_name + "' not found");
      value_map[dep_name] = it->second;
    }

    // Step 2: Process body initializers (constants)
    for (int i = 0; i < body.initializer_size(); ++i) {
      const onnx::TensorProto& tensor = body.initializer(i);
      std::string tensor_name = tensor.name();
      DM dm_const = tensor_to_dm(tensor);
      value_map[tensor_name] = MX(dm_const);
    }

    // Step 3: Create MX symbols for Loop body inputs
    // First input: iteration_num (scalar int)
    // Second input: condition (scalar bool)
    // Remaining inputs: loop-carried dependencies
    std::vector<MX> func_inputs;
    std::vector<std::string> input_names;

    casadi_assert(body.input_size() >= 2,
                 "Loop body must have at least 2 inputs (iter_num, cond)");

    for (int i = 0; i < body.input_size(); ++i) {
      const onnx::ValueInfoProto& input = body.input(i);
      std::string input_name = input.name();

      // Skip if already in value_map (initializer or outer scope)
      if (value_map.count(input_name)) {
        continue;
      }

      // Extract shape
      const onnx::TensorShapeProto& shape = input.type().tensor_type().shape();
      casadi_int rows = get_dimension(shape, 0);
      casadi_int cols = (shape.dim_size() > 1) ? get_dimension(shape, 1) : 1;

      // Create MX symbol
      MX mx_input = MX::sym(input_name, rows, cols);
      value_map[input_name] = mx_input;
      func_inputs.push_back(mx_input);
      input_names.push_back(input_name);
    }

    // Add outer scope dependencies as function inputs
    for (const auto& dep_name : outer_deps) {
      func_inputs.push_back(value_map[dep_name]);
      input_names.push_back(dep_name);
    }

    // Step 4: Process body nodes (same as translate_subgraph_to_function)
    for (int i = 0; i < body.node_size(); ++i) {
      const onnx::NodeProto& node = body.node(i);
      std::string op_type = node.op_type();

      // Gather input tensors
      std::vector<MX> node_inputs;
      for (int j = 0; j < node.input_size(); ++j) {
        std::string input_name = node.input(j);
        if (input_name.empty()) {
          node_inputs.push_back(MX());
          continue;
        }
        casadi_assert(value_map.count(input_name),
                     "Unknown input tensor '" + input_name +
                     "' in Loop body '" + function_name + "' node " +
                     std::to_string(i) + " (op_type: " + op_type + ")");
        node_inputs.push_back(value_map[input_name]);
      }

      // Compute output using shared operation logic
      MX output = process_node_operation(op_type, node, node_inputs);

      // Store output
      casadi_assert(node.output_size() >= 1,
                   "Node must have at least one output");
      std::string output_name = node.output(0);
      value_map[output_name] = output;
    }

    // Step 5: Extract body outputs
    // Outputs: [cond_out, loop_carried_deps..., scan_outputs...]
    std::vector<MX> func_outputs;
    std::vector<std::string> output_names;

    for (int i = 0; i < body.output_size(); ++i) {
      const onnx::ValueInfoProto& output_info = body.output(i);
      std::string output_name = output_info.name();

      casadi_assert(value_map.count(output_name),
                   "Loop body output '" + output_name + "' not found in value_map");

      func_outputs.push_back(value_map[output_name]);
      output_names.push_back(output_name);
    }

    // Step 6: Create and return Function
    return Function(function_name, func_inputs, func_outputs,
                   input_names, output_names);
  }

} // namespace casadi
/// \endcond
