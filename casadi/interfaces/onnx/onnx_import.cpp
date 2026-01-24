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

  // Helper functions for graph processing
  // Process ONNX initializers (pre-loaded constants)
  void process_graph_initializers(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      const OnnxTranslator& translator,
      bool verbose) {

    for (int i = 0; i < graph.initializer_size(); ++i) {
      const onnx::TensorProto& tensor = graph.initializer(i);
      std::string tensor_name = tensor.name();

      if (verbose) {
        uout() << "  Processing initializer: " << tensor_name << std::endl;
      }

      // Convert TensorProto to DM, then to MX
      DM dm_const = translator.tensor_to_dm(tensor);
      value_map[tensor_name] = MX(dm_const);
    }
  }

  // Create MX symbols for graph inputs
  void process_graph_inputs(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      std::vector<MX>& func_inputs,
      std::vector<std::string>& input_names,
      const OnnxTranslator& translator,
      bool verbose) {

    for (int i = 0; i < graph.input_size(); ++i) {
      const onnx::ValueInfoProto& input = graph.input(i);
      std::string input_name = input.name();

      // Skip if already in value_map (it's an initializer, not a variable)
      if (value_map.count(input_name)) {
        if (verbose) {
          uout() << "  Skipping input '" << input_name
                 << "' (it's an initializer)" << std::endl;
        }
        continue;
      }

      // Extract shape
      const onnx::TensorShapeProto& shape =
          input.type().tensor_type().shape();
      casadi_int rows = translator.get_dimension(shape, 0);
      casadi_int cols = translator.get_dimension(shape, 1);

      if (verbose) {
        uout() << "  Creating input: " << input_name
               << " [" << rows << ", " << cols << "]" << std::endl;
      }

      // Create MX symbol
      MX mx_input = MX::sym(input_name, rows, cols);
      value_map[input_name] = mx_input;
      func_inputs.push_back(mx_input);
      input_names.push_back(input_name);
    }
  }

  // Process all nodes in the graph
  void process_graph_nodes(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      OnnxTranslator& translator,
      bool verbose) {

    for (int i = 0; i < graph.node_size(); ++i) {
      const onnx::NodeProto& node = graph.node(i);
      std::string op_type = node.op_type();

      if (verbose) {
        uout() << "  Processing node " << i << ": " << op_type << std::endl;
      }

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
                      "' required by node " + std::to_string(i) +
                      " (op_type: " + op_type + ")");
        node_inputs.push_back(value_map[input_name]);
      }

      // Compute output based on operation type
      MX output;

      // ========== Special Handling: Multi-output operations ==========
      if (op_type == "Split") {
        casadi_assert(node_inputs.size() >= 1, "Split requires 1 input");

        // Get axis attribute (default 0)
        int axis = 0;
        for (int a = 0; a < node.attribute_size(); ++a) {
          if (node.attribute(a).name() == "axis") {
            axis = node.attribute(a).i();
            break;
          }
        }

        // Get split sizes if specified
        std::vector<casadi_int> split_sizes;
        for (int a = 0; a < node.attribute_size(); ++a) {
          if (node.attribute(a).name() == "split") {
            for (int k = 0; k < node.attribute(a).ints_size(); ++k) {
              split_sizes.push_back(node.attribute(a).ints(k));
            }
            break;
          }
        }

        // Perform split based on axis
        std::vector<MX> outputs;
        if (axis == 0) {
          if (split_sizes.empty()) {
            // Equal split - determine from number of outputs
            casadi_int n_outputs = node.output_size();
            casadi_int total_size = node_inputs[0].size1();
            casadi_int size_per = total_size / n_outputs;

            // Build offset vector: [0, size_per, 2*size_per, ..., total_size]
            std::vector<casadi_int> offset;
            offset.push_back(0);
            for (casadi_int j = 0; j < n_outputs; ++j) {
              offset.push_back(offset.back() + size_per);
            }
            outputs = vertsplit(node_inputs[0], offset);
          } else {
            // Build offset vector from split sizes
            std::vector<casadi_int> offset;
            offset.push_back(0);
            for (casadi_int sz : split_sizes) {
              offset.push_back(offset.back() + sz);
            }
            outputs = vertsplit(node_inputs[0], offset);
          }
        } else if (axis == 1) {
          if (split_sizes.empty()) {
            // Equal split - determine from number of outputs
            casadi_int n_outputs = node.output_size();
            casadi_int total_size = node_inputs[0].size2();
            casadi_int size_per = total_size / n_outputs;

            // Build offset vector: [0, size_per, 2*size_per, ..., total_size]
            std::vector<casadi_int> offset;
            offset.push_back(0);
            for (casadi_int j = 0; j < n_outputs; ++j) {
              offset.push_back(offset.back() + size_per);
            }
            outputs = horzsplit(node_inputs[0], offset);
          } else {
            // Build offset vector from split sizes
            std::vector<casadi_int> offset;
            offset.push_back(0);
            for (casadi_int sz : split_sizes) {
              offset.push_back(offset.back() + sz);
            }
            outputs = horzsplit(node_inputs[0], offset);
          }
        } else {
          casadi_error("Split: only axis 0 and 1 supported");
        }

        // Store all outputs
        for (casadi_int j = 0; j < outputs.size(); ++j) {
          value_map[node.output(j)] = outputs[j];
        }
        continue;  // Don't use standard output handling

      // ========== Control Flow Operations (not supported) ==========
      } else if (op_type == "If") {
        casadi_error("ONNX import: 'If' control flow operator is not supported.");
      } else if (op_type == "Loop") {
        casadi_error("ONNX import: 'Loop' control flow operator is not supported.");
      } else if (op_type == "Scan") {
        casadi_error("ONNX import: 'Scan' control flow operator is not supported.");

      // ========== Standard Operations (delegated to helper) ==========
      } else {
        // Check if this is a function call (node has non-empty domain)
        std::string node_domain = node.domain();
        if (!node_domain.empty()) {
          // This is a function call - look up the function definition
          const onnx::ModelProto& model = translator.model_;
          const onnx::FunctionProto* func_proto = nullptr;

          for (int f = 0; f < model.functions_size(); ++f) {
            if (model.functions(f).name() == op_type &&
                model.functions(f).domain() == node_domain) {
              func_proto = &model.functions(f);
              break;
            }
          }

          if (func_proto != nullptr) {
            if (verbose) {
              uout() << "    Function call to: " << node_domain << "." << op_type << std::endl;
            }

            // Create a temporary GraphProto from the FunctionProto
            // This is a bit of a hack but allows us to reuse the graph processing code
            onnx::GraphProto func_graph;
            func_graph.set_name(func_proto->name());

            // Copy nodes from function to graph
            for (int n = 0; n < func_proto->node_size(); ++n) {
              *func_graph.add_node() = func_proto->node(n);
            }

            // Create graph inputs from function inputs
            for (int n = 0; n < func_proto->input_size(); ++n) {
              onnx::ValueInfoProto* input = func_graph.add_input();
              input->set_name(func_proto->input(n));
              // Set default type (double scalar) - actual shapes come from caller
              onnx::TypeProto* type = input->mutable_type();
              onnx::TypeProto::Tensor* tensor_type = type->mutable_tensor_type();
              tensor_type->set_elem_type(onnx::TensorProto::DOUBLE);
            }

            // Create graph outputs from function outputs
            for (int n = 0; n < func_proto->output_size(); ++n) {
              onnx::ValueInfoProto* output_info = func_graph.add_output();
              output_info->set_name(func_proto->output(n));
              onnx::TypeProto* type = output_info->mutable_type();
              onnx::TypeProto::Tensor* tensor_type = type->mutable_tensor_type();
              tensor_type->set_elem_type(onnx::TensorProto::DOUBLE);
            }

            // Build a local value_map for the function scope
            // Map function inputs to the actual values from the call site
            std::map<std::string, MX> func_value_map;
            for (size_t n = 0; n < node_inputs.size() && n < func_proto->input_size(); ++n) {
              func_value_map[func_proto->input(n)] = node_inputs[n];
            }

            // Process function nodes
            process_graph_nodes(func_graph, func_value_map, translator, false);

            // Get function outputs and store them in value_map
            for (int n = 0; n < node.output_size() && n < func_proto->output_size(); ++n) {
              std::string func_output_name = func_proto->output(n);
              casadi_assert(func_value_map.count(func_output_name),
                           "Function output '" + func_output_name + "' not found");
              value_map[node.output(n)] = func_value_map[func_output_name];

              if (verbose) {
                uout() << "      -> " << node.output(n) << std::endl;
              }
            }

            // Skip the normal single-output handling
            continue;
          }
          // Fall through to error if function not found
        }

        // All other operations handled by process_node_operation helper
        // This includes: Add, Sub, Mul, Div, Pow, Sin, Cos, Tan, Asin, Acos, Atan,
        // Sinh, Cosh, Tanh, Asinh, Acosh, Atanh, Exp, Log, Sqrt, Ceil, Floor,
        // Abs, Sign, Neg, Erf, Identity, MatMul, Constant
        output = translator.process_node_operation(op_type, node, node_inputs);
      }

      // Store output (assume single output for now)
      casadi_assert(node.output_size() >= 1,
                    "Node must have at least one output");

      std::string output_name = node.output(0);
      value_map[output_name] = output;

      if (verbose) {
        uout() << "    -> " << output_name << std::endl;
      }
    }
  }

  // Collect graph outputs from value_map
  void collect_graph_outputs(
      const onnx::GraphProto& graph,
      const std::map<std::string, MX>& value_map,
      std::vector<MX>& func_outputs,
      std::vector<std::string>& output_names,
      bool verbose) {

    for (int i = 0; i < graph.output_size(); ++i) {
      const onnx::ValueInfoProto& output = graph.output(i);
      std::string output_name = output.name();

      casadi_assert(value_map.count(output_name),
                    "Unknown output tensor: " + output_name +
                    ". This usually means the ONNX graph contains unsupported operations.");

      func_outputs.push_back(value_map.at(output_name));
      output_names.push_back(output_name);

      if (verbose) {
        uout() << "  Graph output: " << output_name << std::endl;
      }
    }
  }

  Function OnnxTranslator::create(const std::string& name) {
    casadi_assert(has_model_, "No ONNX model loaded. Call load() first.");

    const onnx::GraphProto& graph = model_.graph();

    // Initialize data structures
    std::map<std::string, MX> value_map;
    std::vector<MX> inputs;
    std::vector<std::string> input_names;
    std::vector<MX> outputs;
    std::vector<std::string> output_names;

    if (verbose_) {
    uout() << "Creating CasADi Function from ONNX model: "
           << graph.name() << std::endl;
    uout() << "  Graph has " << graph.initializer_size() << " initializers, "
           << graph.input_size() << " inputs, "
           << graph.output_size() << " outputs, "
           << graph.node_size() << " nodes" << std::endl;
    }

    // Process graph using shared helpers
    process_graph_initializers(graph, value_map, *this, verbose_);
    process_graph_inputs(graph, value_map, inputs, input_names, *this, verbose_);
    process_graph_nodes(graph, value_map, *this, verbose_);
    collect_graph_outputs(graph, value_map, outputs, output_names, verbose_);

    // Create and return CasADi Function
    Function f(name, inputs, outputs, input_names, output_names);

    if (verbose_) {
    uout() << "Created CasADi Function: " << name << std::endl;
    uout() << "  Inputs: " << f.n_in() << std::endl;
    uout() << "  Outputs: " << f.n_out() << std::endl;
    }

    return f;
  }

} // namespace casadi
/// \endcond
