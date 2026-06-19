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

  // ONNX initializers are pre-loaded constants
  void Onnx::process_graph_initializers(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      bool verbose) const {

    for (int i = 0; i < graph.initializer_size(); ++i) {
      const onnx::TensorProto& tensor = graph.initializer(i);
      std::string tensor_name = tensor.name();

      if (verbose) {
        uout() << "  Processing initializer: " << tensor_name << std::endl;
      }

      value_map[tensor_name] = MX(tensor_to_dm(tensor));
    }
  }

  // Create MX symbols for graph inputs
  void Onnx::process_graph_inputs(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      std::vector<MX>& func_inputs,
      std::vector<std::string>& input_names,
      bool verbose) const {

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
      casadi_int rows = get_dimension(shape, 0);
      casadi_int cols = get_dimension(shape, 1);

      if (verbose) {
        uout() << "  Creating input: " << input_name
               << " [" << rows << ", " << cols << "]" << std::endl;
      }

      // Sparsity overlay: if a sparse_initializer shares this input's name it carries a non-dense
      // input pattern -- make a SPARSE input symbol and feed it AS-IS, so CasADi propagates the
      // true sparsity through the body (no densification on entry; densify only where an ONNX op
      // requires
      // the dense layout, e.g. Reshape). Transpose-rep: the ONNX input is declared (c,r).
      Sparsity ov = input_pattern(graph, input_name);
      MX mx_input = ov.is_null() ? MX::sym(input_name, cols, rows) : MX::sym(input_name, ov);
      value_map[input_name] = mx_input;
      func_inputs.push_back(mx_input);
      input_names.push_back(input_name);
    }
  }

  // Process all nodes in the graph
  void Onnx::process_graph_nodes(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      bool verbose) {

    // Per-graph capture of kron operands, keyed by the exporter's node-name tag "kron<k>".
    std::map<std::string, MX> kron_operand_a_, kron_operand_b_;

    for (int i = 0; i < graph.node_size(); ++i) {
      const onnx::NodeProto& node = graph.node(i);
      std::string op_type = node.op_type();

      if (verbose) {
        uout() << "  Processing node " << i << ": " << op_type << std::endl;
      }

      // ========== Special Handling: Kronecker product group ==========
      // The exporter tags its kron node group with node name "kron<k>" (see OP_KRON in
      // onnx_operations.cpp). The 4-D intermediate tensors cannot be MX (MX is 2-D), so we do NOT
      // process the group generically (the helper tensors kron<k>_A4/_B4/_P are never added to the
      // value_map). We capture the two logical operands at the leading Reshapes (whose single data
      // input is the imported, transposed-back A / B) and reconstruct kron(A_imp, B_imp) at the
      // final Reshape. CasADi kron propagates sparsity natively, so the pattern round-trips for
      // free
      // (kron commutes with transpose). This is handled BEFORE the generic input gathering because
      // the helper tensors are intentionally absent from value_map.
      if (node.name().rfind("kron", 0) == 0) {
        const std::string& nm = node.name();  // "kron<k>_<role>"
        std::size_t us = nm.rfind('_');
        std::string tag = nm.substr(0, us), role = nm.substr(us + 1);
        if (role == "A4") {
          kron_operand_a_[tag] = value_map.at(node.input(0));  // data input = imported A
          continue;
        }
        if (role == "B4") {
          kron_operand_b_[tag] = value_map.at(node.input(0));  // data input = imported B
          continue;
        }
        if (role == "P") continue;  // intermediate broadcast product; nothing to record
        if (role == "R") {
          // Final Reshape: emit kron of the captured operands and skip generic 4-D handling.
          casadi_assert(kron_operand_a_.count(tag) && kron_operand_b_.count(tag),
                        "ONNX import: incomplete kron node group '" + tag + "'");
          value_map[node.output(0)] = MX::kron(kron_operand_a_.at(tag), kron_operand_b_.at(tag));
          continue;
        }
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
        casadi_int axis = get_int_attribute(node, "axis", 0);
        casadi_assert(axis == 0 || axis == 1, "Split: only axis 0 and 1 supported");

        // Split sizes: 2nd input (opset 13+), else 'split' attribute (opset <13), else equal
        std::vector<casadi_int> split_sizes;
        if (node_inputs.size() >= 2) {
          DM split_dm = DM(node_inputs[1]);
          for (casadi_int k = 0; k < split_dm.nnz(); ++k) {
            split_sizes.push_back(static_cast<casadi_int>(static_cast<double>(split_dm.nz(k))));
          }
        } else {
          for (int a = 0; a < node.attribute_size(); ++a) {
            if (node.attribute(a).name() == "split") {
              for (int k = 0; k < node.attribute(a).ints_size(); ++k) {
                split_sizes.push_back(node.attribute(a).ints(k));
              }
              break;
            }
          }
        }
        if (split_sizes.empty()) {  // equal split across the outputs
          // Transpose-rep: ONNX axis 0 is CasADi columns, axis 1 is CasADi rows
          casadi_int total = (axis == 0) ? node_inputs[0].size2() : node_inputs[0].size1();
          split_sizes.assign(node.output_size(), total / node.output_size());
        }

        // Offsets [0, s0, s0+s1, ...]; ONNX axis 0 -> CasADi horzsplit, axis 1 -> vertsplit
        std::vector<casadi_int> offset = {0};
        for (casadi_int sz : split_sizes) offset.push_back(offset.back() + sz);
        std::vector<MX> outputs = (axis == 0) ? horzsplit(node_inputs[0], offset)
                                              : vertsplit(node_inputs[0], offset);

        for (casadi_int j = 0; j < outputs.size(); ++j) value_map[node.output(j)] = outputs[j];
        continue;  // Don't use standard output handling

      // ========== Scan: iterate a body over columns -> CasADi Map ==========
      } else if (op_type == "Scan") {
        const onnx::GraphProto* body = get_graph_attribute(node, "body");
        casadi_assert(body != nullptr, "Scan node requires a 'body' subgraph");
        casadi_int num_scan_inputs = get_int_attribute(node, "num_scan_inputs", node.input_size());
        casadi_int M = node.input_size() - num_scan_inputs;  // state variables (reduce_out count)

        // A body that references outer-scope names (captures) signals reduce_in
        std::set<std::string> body_in_names, defined;
        for (int b = 0; b < body->input_size(); ++b) {
          body_in_names.insert(body->input(b).name());
          defined.insert(body->input(b).name());
        }
        for (int nn = 0; nn < body->node_size(); ++nn)
          for (int oo = 0; oo < body->node(nn).output_size(); ++oo)
            defined.insert(body->node(nn).output(oo));
        bool has_capture = false;
        for (int nn = 0; nn < body->node_size() && !has_capture; ++nn)
          for (int ii = 0; ii < body->node(nn).input_size(); ++ii) {
            const std::string& in = body->node(nn).input(ii);
            if (!in.empty() && !defined.count(in)) { has_capture = true; break; }
          }

        if (M == 0 && !has_capture) {
          // Plain map: the 3-D lift reshapes pass through on import (CasADi is 2-D), so
          // node_inputs are the original (rows, n*c) map inputs. Rebuild base.map(n).
          Function base = function_from_graph(*body, op_type + "_body");
          casadi_int n = node_inputs[0].size2() / base.size2_in(0);
          std::vector<MX> outputs =
              base.map(n)(std::vector<MX>(node_inputs.begin(), node_inputs.end()));
          for (int j = 0; j < node.output_size(); ++j) value_map[node.output(j)] = outputs[j];
          continue;
        }

        // Reduce-map: the body wraps a single base call; recover base + masks and rebuild
        // base.map(n, reduce_in, reduce_out).
        const onnx::NodeProto* call = nullptr;
        for (int nn = 0; nn < body->node_size(); ++nn)
          if (!body->node(nn).domain().empty()) { call = &body->node(nn); break; }
        casadi_assert(call, "reduce-map Scan body must contain a base function call");
        casadi_int nin = call->input_size(), nout = call->output_size();

        // reduce_in: base arg is a capture; reduce_out: base result feeds an accumulator Add
        std::vector<bool> reduce_in(nin), reduce_out(nout, false);
        for (casadi_int j = 0; j < nin; ++j) reduce_in[j] = !body_in_names.count(call->input(j));
        for (int nn = 0; nn < body->node_size(); ++nn) {
          if (body->node(nn).op_type() != "Add") continue;
          for (int ii = 0; ii < body->node(nn).input_size(); ++ii)
            for (casadi_int j = 0; j < nout; ++j)
              if (body->node(nn).input(ii) == call->output(j)) reduce_out[j] = true;
        }

        // Base input shapes: captures from the outer tensor, scanned ones from the body input
        std::vector<std::pair<casadi_int, casadi_int>> in_shapes(nin);
        for (casadi_int j = 0; j < nin; ++j) {
          if (reduce_in[j]) {
            MX cap = value_map.at(call->input(j));
            in_shapes[j] = {cap.size1(), cap.size2()};
          } else {
            for (int b = 0; b < body->input_size(); ++b)
              if (body->input(b).name() == call->input(j)) {
                const onnx::TensorShapeProto& sh = body->input(b).type().tensor_type().shape();
                in_shapes[j] = {get_dimension(sh, 0), get_dimension(sh, 1)};
                break;
              }
          }
        }

        const onnx::FunctionProto* fp = find_function(call->op_type(), call->domain());
        casadi_assert(fp, "reduce-map base function '" + call->op_type() + "' not found");
        Function base = function_from_function_proto(*fp, in_shapes, op_type + "_base");

        // base.map args (base order): captures broadcast, repeated come from scan node inputs
        std::vector<MX> args(nin);
        casadi_int n = 0, r = 0;
        for (casadi_int j = 0; j < nin; ++j) {
          if (reduce_in[j]) {
            args[j] = value_map.at(call->input(j));
          } else {
            args[j] = node_inputs[M + r];
            if (n == 0) n = args[j].size2() / base.size2_in(j);
            ++r;
          }
        }
        casadi_assert(n > 0, "reduce-map import: could not infer map size");
        std::vector<MX> outs = base.map(n, reduce_in, reduce_out)(args);

        // Node outputs are [state accumulators (reduce_out), scan outputs (repeated)], base order
        casadi_int si = 0, ci = 0;
        for (casadi_int j = 0; j < nout; ++j) {
          if (reduce_out[j]) value_map[node.output(si++)] = outs[j];
          else
            value_map[node.output(M + ci++)] = outs[j];
        }
        continue;

      // ========== If: two captured branches combined with if_else ==========
      } else if (op_type == "If") {
        const onnx::GraphProto* then_b = get_graph_attribute(node, "then_branch");
        const onnx::GraphProto* else_b = get_graph_attribute(node, "else_branch");
        casadi_assert(then_b && else_b, "If node requires then_branch and else_branch");
        std::vector<MX> t = eval_captured_subgraph(*then_b, value_map);
        std::vector<MX> e = eval_captured_subgraph(*else_b, value_map);
        MX cond = node_inputs[0];
        for (int j = 0; j < node.output_size(); ++j) {
          value_map[node.output(j)] = if_else(cond, t[j], e[j]);
        }
        continue;

      // ========== Control Flow Operations (not supported) ==========
      } else if (op_type == "Loop") {
        casadi_error("ONNX import: 'Loop' control flow operator is not supported.");

      // ========== Standard Operations (delegated to helper) ==========
      } else {
        // A non-empty domain marks a function call.
        std::string node_domain = node.domain();
        if (!node_domain.empty()) {
          const onnx::FunctionProto* func_proto = find_function(op_type, node_domain);

          if (func_proto != nullptr) {
            if (verbose) {
              uout() << "    Function call to: " << node_domain << "." << op_type << std::endl;
            }

            // Inline the function body: only its nodes are needed (process_graph_nodes reads
            // graph.node() only); formal inputs/outputs are wired through func_value_map below.
            onnx::GraphProto func_graph;
            func_graph.set_name(func_proto->name());
            for (int n = 0; n < func_proto->node_size(); ++n) {
              *func_graph.add_node() = func_proto->node(n);
            }

            // Map the function's formal inputs to the call-site values.
            std::map<std::string, MX> func_value_map;
            for (size_t n = 0; n < node_inputs.size() && n < func_proto->input_size(); ++n) {
              func_value_map[func_proto->input(n)] = node_inputs[n];
            }

            process_graph_nodes(func_graph, func_value_map, false);

            // Store the function outputs in the outer value_map.
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

        // All other (single-output) operations
        output = process_node_operation(op_type, node, node_inputs);
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
  void Onnx::collect_graph_outputs(
      const onnx::GraphProto& graph,
      const std::map<std::string, MX>& value_map,
      std::vector<MX>& func_outputs,
      std::vector<std::string>& output_names,
      bool verbose) const {

    for (int i = 0; i < graph.output_size(); ++i) {
      const onnx::ValueInfoProto& output = graph.output(i);
      std::string output_name = output.name();

      casadi_assert(value_map.count(output_name),
                    "Unknown output tensor: " + output_name +
                    ". This usually means the ONNX graph contains unsupported operations.");

      // No output overlay: each op restores its own output sparsity from seeds during export, so
      // the value already carries the exact CasADi pattern (re-propagated natively from the seeds).
      func_outputs.push_back(value_map.at(output_name));
      output_names.push_back(output_name);

      if (verbose) {
        uout() << "  Graph output: " << output_name << std::endl;
      }
    }
  }

  const onnx::FunctionProto* Onnx::find_function(const std::string& name,
                                                 const std::string& domain) const {
    for (int f = 0; f < model_.functions_size(); ++f)
      if (model_.functions(f).name() == name && model_.functions(f).domain() == domain)
        return &model_.functions(f);
    return nullptr;
  }

  Function Onnx::create(const std::string& name) {
    casadi_assert(has_model_, "No ONNX model loaded. Call load() first.");
    return function_from_graph(model_.graph(), name);
  }

  Function Onnx::function_from_graph(const onnx::GraphProto& graph,
                                              const std::string& name) {
    std::map<std::string, MX> value_map;
    std::vector<MX> inputs, outputs;
    std::vector<std::string> input_names, output_names;

    if (verbose_) {
      uout() << "Building CasADi Function '" << name << "' from graph '" << graph.name()
             << "' (" << graph.initializer_size() << " initializers, "
             << graph.input_size() << " inputs, " << graph.output_size() << " outputs, "
             << graph.node_size() << " nodes)" << std::endl;
    }

    process_graph_initializers(graph, value_map, verbose_);
    process_graph_inputs(graph, value_map, inputs, input_names, verbose_);
    process_graph_nodes(graph, value_map, verbose_);
    collect_graph_outputs(graph, value_map, outputs, output_names, verbose_);

    return Function(name, inputs, outputs, input_names, output_names);
  }

  Function Onnx::function_from_function_proto(
      const onnx::FunctionProto& fp,
      const std::vector<std::pair<casadi_int, casadi_int>>& in_shapes,
      const std::string& name) {
    // A FunctionProto carries no shapes; wrap it as a graph with caller-supplied input shapes
    onnx::GraphProto g;
    g.set_name(fp.name());
    for (int j = 0; j < fp.input_size(); ++j) {
      onnx::ValueInfoProto* vi = g.add_input();
      vi->set_name(fp.input(j));
      onnx::TypeProto::Tensor* tt = vi->mutable_type()->mutable_tensor_type();
      tt->set_elem_type(real_type());
      tt->mutable_shape()->add_dim()->set_dim_value(in_shapes[j].first);
      tt->mutable_shape()->add_dim()->set_dim_value(in_shapes[j].second);
    }
    for (int k = 0; k < fp.node_size(); ++k) *g.add_node() = fp.node(k);
    for (int j = 0; j < fp.output_size(); ++j) g.add_output()->set_name(fp.output(j));
    return function_from_graph(g, name);
  }

  std::vector<MX> Onnx::eval_captured_subgraph(const onnx::GraphProto& graph,
                                                        std::map<std::string, MX> scope) {
    // No formal inputs: the branch captures outer tensors, already present in `scope`
    process_graph_initializers(graph, scope, verbose_);
    process_graph_nodes(graph, scope, verbose_);
    std::vector<MX> outputs;
    for (int i = 0; i < graph.output_size(); ++i)
      outputs.push_back(scope.at(graph.output(i).name()));
    return outputs;
  }


  // Convert a constant tensor input to a vector of integers
  static std::vector<casadi_int> constant_ints(const MX& m) {
    casadi_assert(m.is_constant(), "Expected a constant integer tensor");
    DM dm = static_cast<DM>(m);
    std::vector<casadi_int> v;
    for (casadi_int k = 0; k < dm.numel(); ++k)
      v.push_back(static_cast<casadi_int>(dm(k).scalar()));
    return v;
  }

  MX Onnx::process_node_operation(
      const std::string& op_type,
      const onnx::NodeProto& node,
      const std::vector<MX>& node_inputs) {

    MX output;

    // Simple ops via the centralized lookup table
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
          default: break;
        }
      }
    }

    // Everything else needs attributes, multiple inputs, etc.

    if (op_type == "Less") {
      casadi_assert(node_inputs.size() >= 2, "Less requires 2 inputs");
      output = if_else(node_inputs[0] < node_inputs[1], MX(1.0), MX(0.0));

    } else if (op_type == "Equal") {
      casadi_assert(node_inputs.size() >= 2, "Equal requires 2 inputs");
      output = !ne(node_inputs[0], node_inputs[1]);  // CasADi has ne() but no eq()

    } else if (op_type == "LessOrEqual") {
      casadi_assert(node_inputs.size() >= 2, "LessOrEqual requires 2 inputs");
      output = if_else(node_inputs[0] <= node_inputs[1], MX(1.0), MX(0.0));

    } else if (op_type == "Min") {
      casadi_assert(node_inputs.size() >= 2, "Min requires 2 inputs");
      output = fmin(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Max") {
      casadi_assert(node_inputs.size() >= 2, "Max requires 2 inputs");
      output = fmax(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Mod") {
      casadi_assert(node_inputs.size() >= 2, "Mod requires 2 inputs");
      output = fmod(node_inputs[0], node_inputs[1]);

    } else if (op_type == "ReduceSum") {
      // Sum all elements (ReduceMin/Max/L1/L2 are in the op_map table)
      casadi_assert(node_inputs.size() >= 1, "ReduceSum requires 1 input");
      output = sum1(sum2(node_inputs[0]));

    } else if (op_type == "Not") {
      casadi_assert(node_inputs.size() >= 1, "Not requires 1 input");
      output = logic_not(node_inputs[0]);

    } else if (op_type == "And") {
      casadi_assert(node_inputs.size() >= 2, "And requires 2 inputs");
      output = logic_and(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Or") {
      casadi_assert(node_inputs.size() >= 2, "Or requires 2 inputs");
      output = logic_or(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Where") {
      casadi_assert(node_inputs.size() >= 3, "Where requires 3 inputs");
      output = if_else(node_inputs[0], node_inputs[1], node_inputs[2]);

    } else if (op_type == "Identity") {
      casadi_assert(node_inputs.size() >= 1, "Identity requires 1 input");
      output = node_inputs[0];

    } else if (op_type == "Cast") {
      // Everything is double in CasADi (converted at the tensor_to_dm boundary), so Cast is an
      // identity during symbolic computation; the 'to' attribute is irrelevant here.
      casadi_assert(node_inputs.size() >= 1, "Cast requires 1 input");
      output = node_inputs[0];

    } else if (op_type == "MatMul") {
      // Transpose-rep: ONNX MatMul(P,Q)=P*Q is the transpose-rep of q*p, so swap operands.
      casadi_assert(node_inputs.size() >= 2, "MatMul requires 2 inputs");
      output = mtimes(node_inputs[1], node_inputs[0]);

    } else if (op_type == "Gemm") {
      // Transpose-rep: the stored result is (opQ(Q)^T)*(opP(P)^T); recover the CasADi value by
      // swapping operands and the transA/transB roles.
      casadi_assert(node_inputs.size() >= 2, "Gemm requires at least 2 inputs");
      MX A = get_int_attribute(node, "transB", 0) ? node_inputs[1].T() : node_inputs[1];
      MX B = get_int_attribute(node, "transA", 0) ? node_inputs[0].T() : node_inputs[0];
      output = get_float_attribute(node, "alpha", 1.0) * mtimes(A, B);
      if (node_inputs.size() >= 3) {
        output = output + get_float_attribute(node, "beta", 1.0) * node_inputs[2];
      }

    } else if (op_type == "Sum") {
      // Variadic elementwise sum
      casadi_assert(node_inputs.size() >= 1, "Sum requires at least 1 input");
      output = node_inputs[0];
      for (casadi_int idx = 1; idx < node_inputs.size(); ++idx) output = output + node_inputs[idx];

    } else if (op_type == "Pad") {
      // Inverse of the block-diagonal Pad: drop the (dense) block into a zero frame at its offset.
      // Transpose-rep: pads are [axis0_begin, axis1_begin, axis0_end, axis1_end] on the stored
      // (C x R) tensor, i.e. [col_off, row_off, col_end, row_end].
      casadi_assert(node_inputs.size() >= 2, "Pad requires data and pads");
      MX block = node_inputs[0];
      std::vector<casadi_int> pads = constant_ints(node_inputs[1]);
      casadi_int col_off = pads[0], row_off = pads[1];
      casadi_int br = block.size1(), bc = block.size2();
      casadi_int C = col_off + bc + pads[2];
      // Pad with SPARSE zero blocks via horz/vertcat: the block's pattern is preserved and no
      // project/densify is introduced (a triangular block stays triangular), so the block matrix's
      // sparsity emerges natively from CasADi's own propagation.
      MX padded = block;
      if (col_off > 0 || pads[2] > 0) {
        padded = horzcat(MX(Sparsity(br, col_off)), padded, MX(Sparsity(br, pads[2])));
      }
      if (row_off > 0 || pads[3] > 0) {
        padded = vertcat(MX(Sparsity(row_off, C)), padded, MX(Sparsity(pads[3], C)));
      }
      output = padded;

    } else if (op_type == "Einsum") {
      // Reconstruct a CasADi einstein contraction (inverse of the export envelope). The
      // operands arrive as the 2-D reshaped tensors; their column-major vec is the original.
      casadi_assert(node_inputs.size() >= 2, "Einsum requires 2 inputs");
      std::string eq;
      for (char ch : get_string_attribute(node, "equation")) if (ch != ' ') eq += ch;
      size_t comma = eq.find(','), arrow = eq.find("->");
      casadi_assert(comma != std::string::npos && arrow != std::string::npos,
                    "ONNX import: only binary Einsum 'a,b->c' is supported");
      std::string sa = eq.substr(0, comma);
      std::string sb = eq.substr(comma + 1, arrow - comma - 1);
      std::string sc = eq.substr(arrow + 2);

      // Transpose-rep: the exporter reversed each subscript to label the axis-reversed ONNX
      // tensors. node_inputs here are the imported CasADi values (transposed back to LOGICAL axis
      // order), so reverse the subscripts back to natural order to pair them with the logical axes
      // (size1,size2).
      std::reverse(sa.begin(), sa.end());
      std::reverse(sb.begin(), sb.end());
      std::reverse(sc.begin(), sc.end());

      // Each letter's size, read from the operand shapes (axes follow the subscript order). A
      // single-index (vector) operand is stored as a rank-1 ONNX tensor that imports with an
      // ambiguous 1xN/Nx1 orientation, so take its length from numel() (orientation-independent);
      // a 2-index operand reads size1/size2 directly.
      std::map<char, casadi_int> lsize;
      for (int t = 0; t < 2; ++t) {
        const std::string& s = (t == 0) ? sa : sb;
        const MX& m = node_inputs[t];
        if (s.size() == 1) {
          lsize[s[0]] = m.numel();
        } else if (s.size() >= 2) {
          lsize[s[0]] = m.size1();
          lsize[s[1]] = m.size2();
        }
      }
      // Assign a label to each distinct letter
      std::map<char, casadi_int> lab;
      casadi_int next_label = -1;
      for (char ch : sa + sb + sc) if (!lab.count(ch)) lab[ch] = next_label--;

      // Labels in natural order; operands' column-major vec is the original einstein input
      std::vector<casadi_int> da, db, dc, La, Lb, Lc;
      for (char ch : sa) { da.push_back(lsize[ch]); La.push_back(lab[ch]); }
      for (char ch : sb) { db.push_back(lsize[ch]); Lb.push_back(lab[ch]); }
      for (char ch : sc) { dc.push_back(lsize[ch]); Lc.push_back(lab[ch]); }
      output = einstein(vec(node_inputs[0]), vec(node_inputs[1]), da, db, dc, La, Lb, Lc);

    } else if (op_type == "Det") {
      casadi_assert(node_inputs.size() >= 1, "Det requires 1 input");
      output = det(node_inputs[0]);

    } else if (op_type == "ReduceLogSumExp") {
      // log(sum(exp(x))) over all elements
      casadi_assert(node_inputs.size() >= 1, "ReduceLogSumExp requires 1 input");
      output = log(sum1(sum2(exp(node_inputs[0]))));

    } else if (op_type == "Constant") {
      // A dense constant uses the 'value' attribute; a sparse one uses 'sparse_value' (COO)
      const onnx::AttributeProto* value_attr = nullptr;
      const onnx::AttributeProto* sparse_attr = nullptr;
      for (int a = 0; a < node.attribute_size(); ++a) {
        const std::string& an = node.attribute(a).name();
        if (an == "value") value_attr = &node.attribute(a);
        else if (an == "sparse_value") sparse_attr = &node.attribute(a);
      }
      if (sparse_attr != nullptr) {
        output = MX(sparse_tensor_to_dm(sparse_attr->sparse_tensor()));
      } else {
        casadi_assert(value_attr != nullptr,
                     "Constant node must have a 'value' or 'sparse_value' attribute");
        output = MX(tensor_to_dm(value_attr->t()));
      }

    // Complex tensor operations (Transpose is handled by the op_map table)
    } else if (op_type == "Reshape") {
      casadi_assert(node_inputs.size() >= 2,
                    "Reshape operation requires 2 inputs (data and shape)");
      // Second input is the target shape - should be a constant
      casadi_assert(node_inputs[1].is_constant(),
                    "Reshape shape must be a constant");
      DM shape_dm = static_cast<DM>(node_inputs[1]);
      // ONNX Reshape is PSEUDO-DENSE (it operates on the full numel), so densify the operand here:
      // this is the ONE place the imported graph must drop sparsity, because CasADi's own reshape
      // PRESERVES it (a sparse operand would reshape to a flat with only nnz entries, shifting
      // every
      // downstream Gather/Scatter index). True sparsity propagates everywhere else.
      if (shape_dm.numel() > 2) {
        // 3-D reshape (the Map/Scan lift envelope): CasADi is 2-D, so pass through unchanged
        output = densify(node_inputs[0]);
      } else {
        // Transpose-rep: the ONNX target (s0,s1) is the reverse of the CasADi target.
        casadi_int s0 = static_cast<casadi_int>(shape_dm(0).scalar());
        casadi_int s1 = (shape_dm.numel() > 1) ? static_cast<casadi_int>(shape_dm(1).scalar()) : 1;
        output = reshape(densify(node_inputs[0]), s1, s0);
      }

    } else if (op_type == "Concat") {
      // Transpose-rep: ONNX axis 0 is CasADi columns (horzcat), axis 1 is CasADi rows (vertcat)
      casadi_int axis = get_int_attribute(node, "axis", 0);
      if (axis == 0) {
        output = horzcat(node_inputs);
      } else if (axis == 1) {
        output = vertcat(node_inputs);
      } else {
        casadi_error("Concat with axis=" + std::to_string(axis) +
                     " not supported. Only axis=0 (vertcat) and axis=1 (horzcat) are supported.");
      }

    } else if (op_type == "Slice") {
      // Inputs: data, starts, ends, [axes], [steps].
      casadi_assert(node_inputs.size() >= 3,
                    "Slice requires at least 3 inputs (data, starts, ends)");

      MX data = node_inputs[0];
      casadi_assert(node_inputs[1].is_constant() && node_inputs[2].is_constant(),
                    "Slice starts and ends must be constants");

      DM starts_dm = static_cast<DM>(node_inputs[1]);
      DM ends_dm = static_cast<DM>(node_inputs[2]);

      // Axes (default [0, 1, ...]) and steps (default all 1s)
      std::vector<casadi_int> axes, steps;
      if (node_inputs.size() >= 4 && !node_inputs[3].is_empty()) {
        axes = constant_ints(node_inputs[3]);
      } else {
        for (casadi_int k = 0; k < starts_dm.numel(); ++k) axes.push_back(k);
      }
      if (node_inputs.size() >= 5 && !node_inputs[4].is_empty()) {
        steps = constant_ints(node_inputs[4]);
      } else {
        steps.assign(starts_dm.numel(), 1);
      }

      casadi_assert(axes.size() <= 2, "Slice: only up to 2D slicing supported");

      // Transpose-rep: ONNX axis 0 is CasADi columns, axis 1 is CasADi rows. Build a row- and a
      // column-Slice (default = all) from whichever axes are present, with their steps, then index.
      Slice row_slice, col_slice;
      casadi_int nrow = data.size1(), ncol = data.size2();
      for (casadi_int a = 0; a < static_cast<casadi_int>(axes.size()); ++a) {
        casadi_int st = static_cast<casadi_int>(starts_dm(a).scalar());
        casadi_int en = static_cast<casadi_int>(ends_dm(a).scalar());
        casadi_int sp = steps[a];
        if (axes[a] == 0) {
          col_slice = Slice(st, en > ncol ? ncol : en, sp);  // -> CasADi cols
        } else {
          row_slice = Slice(st, en > nrow ? nrow : en, sp);  // -> CasADi rows
        }
      }
      output = data(row_slice, col_slice);

    } else if (op_type == "Gather") {
      // Extract element(s) at the given indices along an axis.
      casadi_assert(node_inputs.size() >= 2, "Gather requires data and indices");

      // Gather indexes the pseudo-dense column-major layout -> densify a sparse operand.
      MX data = densify(node_inputs[0]);
      MX indices_mx = node_inputs[1];
      casadi_int axis = get_int_attribute(node, "axis", 0);

      casadi_assert(indices_mx.is_constant(), "Gather indices must be constant");
      DM indices_dm = static_cast<DM>(indices_mx);

      if (data.size2() == 1) {
        // Transpose-rep getnonzeros: data is a column-flat (numel x 1); gather elements by flat
        // index regardless of the ONNX axis.
        if (indices_dm.numel() == 1) {
          output = data(static_cast<casadi_int>(indices_dm(0).scalar()), Slice());
        } else {
          output = data(constant_ints(indices_mx), Slice());
        }
      } else if (indices_dm.numel() == 1) {
        casadi_int idx = static_cast<casadi_int>(indices_dm(0).scalar());
        if (axis == 0) {
          output = data(idx, Slice());      // a row
        } else if (axis == 1) {
          output = data(Slice(), idx);      // a column
        } else {
          casadi_error("Gather: only axis 0 and 1 supported for 2D tensors");
        }
      } else {
        // Multiple indices: gather the rows (axis 0) / columns (axis 1).
        std::vector<casadi_int> indices = constant_ints(indices_mx);
        if (axis == 0) {
          std::vector<MX> rows;
          for (casadi_int idx : indices) rows.push_back(data(idx, Slice()));
          output = vertcat(rows);
        } else if (axis == 1) {
          std::vector<MX> cols;
          for (casadi_int idx : indices) cols.push_back(data(Slice(), idx));
          output = horzcat(cols);
        } else {
          casadi_error("Gather: only axis 0 and 1 supported for 2D tensors");
        }
      }

    } else if (op_type == "ScatterElements") {
      // data with data.nz[indices] {=,+=} updates -- nz indexing is pseudo-dense, so densify the
      // target. Constant indices -> (set/add)nonzeros; runtime indices -> the _PARAM variants.
      casadi_assert(node_inputs.size() >= 3, "ScatterElements requires data, indices, updates");
      MX se_data = densify(node_inputs[0]);
      MX se_idx = node_inputs[1];
      MX se_upd = vec(node_inputs[2]);
      if (get_string_attribute(node, "reduction") == "add") {
        // {add,addparam}nonzeros: out = data with data.nz[idx] += updates, accumulating DUPLICATE
        // idx. No direct MX builder; it is the reverse-mode adjoint of a (parametric)
        // getnonzeros --
        // jtimes(w.nz[idx], w, updates, true) is the scatter-add (dups summed) into a zero frame.
        MX w = MX::sym("scel_w", se_data.numel(), 1), g;
        if (se_idx.is_constant()) {
          w.get_nz(g, false, Matrix<casadi_int>(constant_ints(se_idx)));
        } else {
          w.get_nz(g, false, vec(se_idx));   // runtime indices -> getnonzeros_param
        }
        output = se_data + reshape(jtimes(g, w, se_upd, true), se_data.size1(), se_data.size2());
      } else {
        output = se_data;
        if (se_idx.is_constant()) {
          output.set_nz(se_upd, false, Matrix<casadi_int>(constant_ints(se_idx)));
        } else {
          output.set_nz(se_upd, false, vec(se_idx));   // runtime indices -> setnonzeros_param
        }
      }

    } else if (op_type == "ScatterND") {
      // setnonzeros: overwrite x's cells at the given (onnx_row, onnx_col) coordinates. With x
      // stored as x^T, coord (a,b) is CasADi cell (b,a) -> column-major position a*nrow+b, which is
      // exactly the original nz position, so reconstruct via set_nz on the 2-D value.
      casadi_assert(node_inputs.size() >= 3, "ScatterND requires data, indices, updates");
      MX data = densify(node_inputs[0]);  // nz indexing is pseudo-dense -> densify the target
      std::vector<casadi_int> coords = constant_ints(node_inputs[1]);  // [r0,c0,r1,c1,...]
      casadi_int R = data.size1();
      std::vector<casadi_int> pos(coords.size() / 2);
      for (casadi_int i = 0; i < static_cast<casadi_int>(pos.size()); ++i) {
        pos[i] = coords[2 * i] * R + coords[2 * i + 1];
      }
      if (get_string_attribute(node, "reduction") == "add") {
        // ADDNONZEROS: out = data with out.nz[pos] += updates, accumulating DUPLICATE pos. An
        // addnonzeros MX cannot be built directly; it is the reverse-mode adjoint of a getnonzeros.
        // Build g = w.nz[pos] on a dense numel-vector w, then jtimes(g, w, updates, /*tr*/true) is
        // exactly the scatter-add of updates into a zero frame (duplicates summed); add it to data.
        MX w = MX::sym("addnz_w", data.numel(), 1), g;
        w.get_nz(g, false, Matrix<casadi_int>(pos));
        MX scatter = jtimes(g, w, vec(node_inputs[2]), true);
        output = data + reshape(scatter, data.size1(), data.size2());
      } else {
        output = data;
        output.set_nz(node_inputs[2], false, Matrix<casadi_int>(pos));
      }

    } else if (op_type == "Tile") {
      // Repeat the tensor along each dimension by repeats = [rows, cols].
      casadi_assert(node_inputs.size() >= 2, "Tile requires data and repeats inputs");
      MX data = node_inputs[0];
      MX repeats_mx = node_inputs[1];
      casadi_assert(repeats_mx.is_constant(), "Tile repeats must be constant");
      DM repeats_dm = static_cast<DM>(repeats_mx);

      casadi_int rows_repeat = 1, cols_repeat = 1;
      if (repeats_dm.numel() >= 1) rows_repeat = static_cast<casadi_int>(repeats_dm(0).scalar());
      if (repeats_dm.numel() >= 2) cols_repeat = static_cast<casadi_int>(repeats_dm(1).scalar());

      // Transpose-rep: ONNX repeats [rows,cols] map to CasADi repmat(data, cols, rows)
      output = repmat(data, cols_repeat, rows_repeat);

    } else if (op_type == "GatherElements") {
      // out[..,j,..] = data[.., idx[..,j,..], ..] along `axis`. The export of OP_GETNONZEROS_PARAM
      // emits this on a 1 x N flat row with a 1 x M index row (-> imported as N x 1 / M x 1 cols),
      // which is exactly the parametric getnonzeros y = data.nz[idx]. Reconstruct via get_nz so the
      // sparse/runtime composition (loc-mapping gather then data gather) rebuilds the right MX.
      casadi_assert(node_inputs.size() >= 2, "GatherElements requires data and indices");
      MX data = densify(node_inputs[0]);  // nz indexing is pseudo-dense -> densify the operand
      MX indices_mx = node_inputs[1];
      // Orient indices as a column vector of nz-indices (a 1 x M ONNX row imports as M x 1).
      MX idx_vec = vec(indices_mx);
      MX flat = vec(data);  // column-major flat of data: flat.nz[k] = data's k-th nonzero
      if (indices_mx.is_constant()) {
        // Constant indices -> the plain constant getnonzeros (OP_GETNONZEROS) path.
        flat.get_nz(output, false, Matrix<casadi_int>(constant_ints(indices_mx)));
      } else {
        // Runtime indices -> parametric getnonzeros (OP_GETNONZEROS_PARAM).
        flat.get_nz(output, false, idx_vec);
      }

    } else {
      casadi_error("Unsupported operation '" + op_type + "'");
    }

    return output;
  }

} // namespace casadi
/// \endcond
