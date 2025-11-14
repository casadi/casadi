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

  Function OnnxTranslator::create(const std::string& name) {
    casadi_assert(has_model_, "No ONNX model loaded. Call load() first.");

    const onnx::GraphProto& graph = model_.graph();

    // Step 1: Initialize data structures
    // Map tensor names to MX expressions (the "symbol table")
    std::map<std::string, MX> value_map;

    // Collect function inputs/outputs
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

    // Step 2: Process graph initializers (pre-loaded constants)
    for (int i = 0; i < graph.initializer_size(); ++i) {
      const onnx::TensorProto& tensor = graph.initializer(i);
      std::string tensor_name = tensor.name();

      if (verbose_) {
        uout() << "  Processing initializer: " << tensor_name << std::endl;
      }

      // Convert TensorProto to DM, then to MX
      DM dm_const = tensor_to_dm(tensor);
      value_map[tensor_name] = MX(dm_const);
    }

    // Step 3: Create MX symbols for graph inputs
    for (int i = 0; i < graph.input_size(); ++i) {
      const onnx::ValueInfoProto& input = graph.input(i);
      std::string input_name = input.name();

      // Skip if already in value_map (it's an initializer, not a variable)
      if (value_map.count(input_name)) {
        if (verbose_) {
          uout() << "  Skipping input '" << input_name
                 << "' (it's an initializer)" << std::endl;
        }
        continue;
      }

      // Extract shape
      const onnx::TensorShapeProto& shape =
          input.type().tensor_type().shape();
      casadi_int rows = get_dimension(shape, 0);
      casadi_int cols = get_dimension(shape, 1);

      if (verbose_) {
        uout() << "  Creating input: " << input_name
               << " [" << rows << ", " << cols << "]" << std::endl;
      }

      // Create MX symbol
      MX mx_input = MX::sym(input_name, rows, cols);
      value_map[input_name] = mx_input;
      inputs.push_back(mx_input);
      input_names.push_back(input_name);
    }

    // Step 4: Process nodes in topological order
    for (int i = 0; i < graph.node_size(); ++i) {
      const onnx::NodeProto& node = graph.node(i);
      std::string op_type = node.op_type();

      if (verbose_) {
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

      // ========== Special Handling: Complex tensor operations ==========
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

        if (verbose_) {
          uout() << "    Slice: basic 2D slicing applied" << std::endl;
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

        // For now, implement a simplified version
        // Full implementation requires advanced CasADi indexing capabilities
        casadi_warning("ONNX import: GatherElements is not fully supported. "
                      "Using first input as placeholder.");

        // Just pass through the data for now
        output = node_inputs[0];

        if (verbose_) {
          uout() << "    GatherElements: simplified implementation (axis=" << axis << ")" << std::endl;
        }

      // ========== Control Flow Operations ==========
      } else if (op_type == "If") {
        // ========== Control Flow: If operator ==========
        // If operator: conditional execution with subgraphs
        casadi_assert(node_inputs.size() >= 1,
                     "If operator requires condition input");

        MX condition = node_inputs[0];

        // Extract then_branch and else_branch from attributes
        const onnx::GraphProto* then_branch = nullptr;
        const onnx::GraphProto* else_branch = nullptr;

        for (int a = 0; a < node.attribute_size(); ++a) {
          const onnx::AttributeProto& attr = node.attribute(a);
          if (attr.name() == "then_branch" && attr.type() == onnx::AttributeProto::GRAPH) {
            then_branch = &attr.g();
          } else if (attr.name() == "else_branch" && attr.type() == onnx::AttributeProto::GRAPH) {
            else_branch = &attr.g();
          }
        }

        casadi_assert(then_branch != nullptr,
                     "If operator must have 'then_branch' attribute");
        casadi_assert(else_branch != nullptr,
                     "If operator must have 'else_branch' attribute");

        if (verbose_) {
          uout() << "    If operator: analyzing branches" << std::endl;
          uout() << "      then_branch: " << then_branch->node_size() << " nodes" << std::endl;
          uout() << "      else_branch: " << else_branch->node_size() << " nodes" << std::endl;
        }

        // Build set of available variables in current scope
        std::set<std::string> available_vars;
        for (const auto& kv : value_map) {
          available_vars.insert(kv.first);
        }

        // Analyze outer scope dependencies for both branches
        std::vector<std::string> then_deps =
            analyze_outer_scope_dependencies(*then_branch, available_vars);
        std::vector<std::string> else_deps =
            analyze_outer_scope_dependencies(*else_branch, available_vars);

        // Take union of dependencies (both branches need access to all deps)
        std::set<std::string> all_deps_set(then_deps.begin(), then_deps.end());
        all_deps_set.insert(else_deps.begin(), else_deps.end());
        std::vector<std::string> all_deps(all_deps_set.begin(), all_deps_set.end());

        if (verbose_ && !all_deps.empty()) {
          uout() << "      Outer scope dependencies: ";
          for (const auto& dep : all_deps) {
            uout() << dep << " ";
          }
          uout() << std::endl;
        }

        // Translate both branches to CasADi Functions
        Function f_then = translate_subgraph_to_function(
            *then_branch, "if_then_" + std::to_string(i), value_map, all_deps);
        Function f_else = translate_subgraph_to_function(
            *else_branch, "if_else_" + std::to_string(i), value_map, all_deps);

        // Verify output compatibility
        casadi_assert(f_then.n_out() == f_else.n_out(),
                     "If branches must have same number of outputs");
        casadi_assert(node.output_size() == f_then.n_out(),
                     "If operator output count must match branch outputs");

        // Create if_else Function using CasADi's Function::if_else
        Function f_if = Function::if_else("if_" + std::to_string(i), f_then, f_else);

        // Prepare inputs: condition + branch inputs + outer scope dependencies
        std::vector<MX> if_all_inputs;

        // First input is the condition
        if_all_inputs.push_back(condition);

        // Add outer scope dependencies as inputs
        for (const auto& dep_name : all_deps) {
          casadi_assert(value_map.count(dep_name),
                       "Outer scope dependency '" + dep_name + "' not found");
          if_all_inputs.push_back(value_map[dep_name]);
        }

        // Call the if_else function with all inputs in one vector
        std::vector<MX> if_outputs = f_if(if_all_inputs);

        // Store outputs
        casadi_assert(if_outputs.size() == node.output_size(),
                     "If operator output size mismatch");

        for (size_t j = 0; j < node.output_size(); ++j) {
          std::string output_name = node.output(j);
          value_map[output_name] = if_outputs[j];
          if (verbose_) {
            uout() << "      -> " << output_name << std::endl;
          }
        }

        // Skip the normal single-output handling below
        continue;

      } else if (op_type == "Loop") {
        // ========== Control Flow: Loop operator ==========
        // Loop operator: iteration with state (like while/for loop)
        casadi_assert(node_inputs.size() >= 2,
                     "Loop requires at least 2 inputs (max_iter, cond)");

        MX max_iter_mx = node_inputs[0];  // May be empty for unbounded
        MX initial_cond = node_inputs[1];

        // Remaining inputs are loop-carried dependencies
        std::vector<MX> loop_carried_initial;
        for (size_t j = 2; j < node_inputs.size(); ++j) {
          loop_carried_initial.push_back(node_inputs[j]);
        }

        // Extract body subgraph
        const onnx::GraphProto* body = nullptr;
        for (int a = 0; a < node.attribute_size(); ++a) {
          const onnx::AttributeProto& attr = node.attribute(a);
          if (attr.name() == "body" && attr.type() == onnx::AttributeProto::GRAPH) {
            body = &attr.g();
            break;
          }
        }

        casadi_assert(body != nullptr, "Loop must have 'body' attribute");

        if (verbose_) {
          uout() << "    Loop operator: analyzing body" << std::endl;
          uout() << "      body: " << body->node_size() << " nodes" << std::endl;
          uout() << "      loop-carried deps: " << loop_carried_initial.size() << std::endl;
        }

        // Build set of available variables in current scope
        std::set<std::string> available_vars;
        for (const auto& kv : value_map) {
          available_vars.insert(kv.first);
        }

        // Analyze outer scope dependencies
        std::vector<std::string> outer_deps =
            analyze_outer_scope_dependencies(*body, available_vars);

        if (verbose_ && !outer_deps.empty()) {
          uout() << "      Outer scope dependencies: ";
          for (const auto& dep : outer_deps) {
            uout() << dep << " ";
          }
          uout() << std::endl;
        }

        // Translate body to Function
        Function body_func = translate_loop_body_to_function(
            *body, "loop_body_" + std::to_string(i), value_map, outer_deps);

        // Determine iteration count
        casadi_int n_iter = 10;  // Default
        if (!max_iter_mx.is_empty()) {
          if (max_iter_mx.is_constant()) {
            DM max_iter_dm = static_cast<DM>(max_iter_mx);
            n_iter = static_cast<casadi_int>(max_iter_dm.scalar());
          } else {
            casadi_warning("Loop: symbolic iteration count not supported, using default 10");
          }
        }

        if (verbose_) {
          uout() << "      Iterations: " << n_iter << std::endl;
        }

        // NOTE: For MVP, we execute fixed iterations
        // TODO: Implement conditional termination using body's cond_out

        // Use mapaccum for stateful iteration
        // Body outputs: [cond_out, loop_carried..., scan_outputs...]
        // We'll use only the loop-carried part for mapaccum state
        Function loop_func = body_func.mapaccum("loop_" + std::to_string(i), n_iter);

        // Prepare inputs: initial condition + loop-carried deps + outer deps
        std::vector<MX> loop_inputs;
        loop_inputs.push_back(initial_cond);  // Initial condition
        loop_inputs.insert(loop_inputs.end(),
                          loop_carried_initial.begin(),
                          loop_carried_initial.end());

        // Add outer scope dependencies
        for (const auto& dep_name : outer_deps) {
          casadi_assert(value_map.count(dep_name),
                       "Outer scope dependency '" + dep_name + "' not found");
          loop_inputs.push_back(value_map[dep_name]);
        }

        // Execute loop
        std::vector<MX> loop_outputs = loop_func(loop_inputs);

        // Store outputs
        // Loop outputs: final loop-carried values + scan outputs
        casadi_assert(loop_outputs.size() >= node.output_size(),
                     "Loop output size mismatch");

        for (size_t j = 0; j < node.output_size(); ++j) {
          std::string output_name = node.output(j);
          value_map[output_name] = loop_outputs[j];
          if (verbose_) {
            uout() << "      -> " << output_name << std::endl;
          }
        }

        // Skip normal single-output handling
        continue;

      } else if (op_type == "Scan") {
        // ========== Control Flow: Scan operator ==========
        // Scan operator: functional iteration with NO outer scope access
        casadi_assert(node_inputs.size() >= 1, "Scan requires at least 1 input");

        // Parse attributes
        casadi_int num_scan_inputs = 0;
        for (const auto& attr : node.attribute()) {
          if (attr.name() == "num_scan_inputs") {
            num_scan_inputs = attr.i();
            break;
          }
        }

        // Split inputs: state variables + scan inputs
        size_t num_state = node_inputs.size() - num_scan_inputs;
        std::vector<MX> initial_state(node_inputs.begin(),
                                       node_inputs.begin() + num_state);
        std::vector<MX> scan_inputs(node_inputs.begin() + num_state,
                                     node_inputs.end());

        // Extract body subgraph
        const onnx::GraphProto* body = nullptr;
        for (const auto& attr : node.attribute()) {
          if (attr.name() == "body" && attr.type() == onnx::AttributeProto::GRAPH) {
            body = &attr.g();
            break;
          }
        }

        casadi_assert(body != nullptr, "Scan must have 'body' attribute");

        if (verbose_) {
          uout() << "    Scan operator: analyzing body" << std::endl;
          uout() << "      body: " << body->node_size() << " nodes" << std::endl;
          uout() << "      state vars: " << initial_state.size() << std::endl;
          uout() << "      scan inputs: " << scan_inputs.size() << std::endl;
        }

        // Verify no outer scope access (Scan restriction)
        std::set<std::string> empty_scope;
        std::vector<std::string> outer_deps =
            analyze_outer_scope_dependencies(*body, empty_scope);

        if (!outer_deps.empty()) {
          std::string deps_str;
          for (const auto& dep : outer_deps) {
            if (!deps_str.empty()) deps_str += ", ";
            deps_str += dep;
          }
          casadi_error("Scan body cannot access outer scope variables. Found: " + deps_str);
        }

        // Translate body to Function
        std::map<std::string, MX> empty_outer_vars;
        Function body_func = translate_subgraph_to_function(
            *body, "scan_body_" + std::to_string(i), empty_outer_vars, {});

        // Determine scan length from first scan input
        casadi_int scan_length = 1;
        if (!scan_inputs.empty()) {
          scan_length = scan_inputs[0].size1();  // Assuming axis=0
        }

        if (verbose_) {
          uout() << "      Scan length: " << scan_length << std::endl;
        }

        // Use Function::map for parallel execution
        Function scan_func = body_func.map(scan_length, "serial");

        // Prepare inputs
        std::vector<MX> scan_func_inputs;
        scan_func_inputs.insert(scan_func_inputs.end(),
                               initial_state.begin(),
                               initial_state.end());
        scan_func_inputs.insert(scan_func_inputs.end(),
                               scan_inputs.begin(),
                               scan_inputs.end());

        // Execute scan
        std::vector<MX> scan_outputs = scan_func(scan_func_inputs);

        // Store outputs
        for (size_t j = 0; j < node.output_size(); ++j) {
          std::string output_name = node.output(j);
          value_map[output_name] = scan_outputs[j];
          if (verbose_) {
            uout() << "      -> " << output_name << std::endl;
          }
        }

        // Skip normal single-output handling
        continue;

      // ========== Standard Operations (delegated to helper) ==========
      } else {
        // All other operations handled by process_node_operation helper
        // This includes: Add, Sub, Mul, Div, Pow, Sin, Cos, Tan, Asin, Acos, Atan,
        // Sinh, Cosh, Tanh, Asinh, Acosh, Atanh, Exp, Log, Sqrt, Ceil, Floor,
        // Abs, Sign, Neg, Erf, Identity, MatMul, Constant
        output = process_node_operation(op_type, node, node_inputs);
      }

      // Store output (assume single output for now)
      casadi_assert(node.output_size() >= 1,
                    "Node must have at least one output");

      std::string output_name = node.output(0);
      value_map[output_name] = output;

      if (verbose_) {
        uout() << "    -> " << output_name << std::endl;
      }
    }

    // Step 5: Collect graph outputs
    for (int i = 0; i < graph.output_size(); ++i) {
      const onnx::ValueInfoProto& output = graph.output(i);
      std::string output_name = output.name();

      casadi_assert(value_map.count(output_name),
                    "Unknown output tensor: " + output_name +
                    ". This usually means the ONNX graph contains unsupported operations.");

      outputs.push_back(value_map[output_name]);
      output_names.push_back(output_name);

      if (verbose_) {
        uout() << "  Graph output: " << output_name << std::endl;
      }
    }

    // Step 6: Create and return CasADi Function
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
