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
#include <casadi/core/casadi_misc.hpp>
#include <casadi/core/calculus.hpp>
#include <casadi/core/mx.hpp>
#include <casadi/config.h>
#include <fstream>

namespace casadi {
  extern "C"
  int CASADI_TRANSLATOR_ONNX_EXPORT
  casadi_register_translator_onnx(TranslatorInternal::Plugin* plugin) {
    plugin->creator = OnnxTranslator::creator;
    plugin->name = "onnx";
    plugin->doc = OnnxTranslator::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &OnnxTranslator::options_;
    return 0;
  }

  extern "C"
  void CASADI_TRANSLATOR_ONNX_EXPORT casadi_load_translator_onnx() {
    TranslatorInternal::registerPlugin(casadi_register_translator_onnx);
  }

  OnnxTranslator::OnnxTranslator()
    : TranslatorInternal("onnx"), has_model_(false) {
  }

  OnnxTranslator::~OnnxTranslator() {
  }

  const Options OnnxTranslator::options_
  = {{&TranslatorInternal::options_},
     {{"opset_version",
       {OT_INT,
        "ONNX opset version to use when exporting (default: 13)"}}}
  };

  void OnnxTranslator::init(const Dict& opts) {
    // Call base class
    TranslatorInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="opset_version") {
        // Store for later use
      }
    }
  }

  void OnnxTranslator::load(const std::string& filename) {
    // Read ONNX file
    std::ifstream input(filename, std::ios::binary);
    casadi_assert(input.good(), "Failed to open ONNX file: " + filename);

    // Parse protobuf
    casadi_assert(model_.ParseFromIstream(&input),
                  "Failed to parse ONNX model from file: " + filename);

    input.close();
    has_model_ = true;

    if (verbose_) {
      uout() << "Loaded ONNX model: " << model_.graph().name() << std::endl;
      uout() << "  Inputs: " << model_.graph().input_size() << std::endl;
      uout() << "  Outputs: " << model_.graph().output_size() << std::endl;
      uout() << "  Nodes: " << model_.graph().node_size() << std::endl;
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

  void OnnxTranslator::set_dimension(const std::string& name, casadi_int dim) {
    dimension_overrides_[name] = dim;

    if (verbose_) {
      uout() << "Set dimension override: " << name << " = " << dim << std::endl;
    }
  }

  casadi_int OnnxTranslator::get_dimension(
      const onnx::TensorShapeProto& shape, int idx) const {
    // Out of range or no dimensions = scalar (1)
    if (idx >= shape.dim_size()) return 1;

    const auto& dim = shape.dim(idx);

    // Concrete dimension
    if (dim.has_dim_value()) {
      return static_cast<casadi_int>(dim.dim_value());
    }

    // Symbolic dimension - look up in overrides
    if (dim.has_dim_param()) {
      std::string param_name = dim.dim_param();

      auto it = dimension_overrides_.find(param_name);
      casadi_assert(it != dimension_overrides_.end(),
                    "Symbolic dimension '" + param_name + "' not specified. " +
                    "Call set_dimension(\"" + param_name + "\", value) before create().");

      return it->second;
    }

    // No dimension info - default to 1
    return 1;
  }

  DM OnnxTranslator::tensor_to_dm(const onnx::TensorProto& tensor) const {
    // Extract shape
    std::vector<casadi_int> dims;
    for (int i = 0; i < tensor.dims_size(); ++i) {
      dims.push_back(static_cast<casadi_int>(tensor.dims(i)));
    }

    if (verbose_) {
      uout() << "      tensor_to_dm: dims_size=" << tensor.dims_size()
             << ", double_data_size=" << tensor.double_data_size()
             << ", raw_data_size=" << tensor.raw_data().size()
             << ", data_type=" << tensor.data_type();
      if (dims.size() > 0) {
        uout() << ", dims[0]=" << dims[0];
      }
      uout() << std::endl;
    }

    // Determine matrix dimensions (assume 2D or less)
    casadi_int rows = dims.size() > 0 ? dims[0] : 1;
    casadi_int cols = dims.size() > 1 ? dims[1] : 1;

    if (verbose_) {
      uout() << "      rows=" << rows << ", cols=" << cols << std::endl;
    }

    // For higher dimensional tensors, flatten to 2D for now
    if (dims.size() > 2) {
      casadi_warning("ONNX tensor has " + std::to_string(dims.size()) +
                     " dimensions, flattening to 2D");
      cols = 1;
      for (size_t i = 1; i < dims.size(); ++i) {
        cols *= dims[i];
      }
    }

    // Create DM
    DM dm(rows, cols);

    // Fill with data (only handle DOUBLE for now)
    casadi_assert(tensor.data_type() == onnx::TensorProto::DOUBLE,
                  "Only DOUBLE tensors supported in import, got type " +
                  std::to_string(tensor.data_type()));

    // ONNX can store data in either double_data field or raw_data field
    // Check raw_data first since that's what our export uses
    if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
      // Data in raw bytes
      const std::string& raw = tensor.raw_data();
      casadi_assert(raw.size() == rows * cols * sizeof(double),
                    "Raw data size mismatch: expected " +
                    std::to_string(rows * cols * sizeof(double)) + ", got " +
                    std::to_string(raw.size()));

      const double* data_ptr = reinterpret_cast<const double*>(raw.data());
      for (casadi_int i = 0; i < rows * cols; ++i) {
        dm(i) = data_ptr[i];
      }
    } else if (tensor.double_data_size() > 0) {
      // Data in typed field
      casadi_assert(tensor.double_data_size() == rows * cols,
                    "Tensor data size mismatch: expected " +
                    std::to_string(rows * cols) + ", got " +
                    std::to_string(tensor.double_data_size()));

      for (casadi_int i = 0; i < tensor.double_data_size(); ++i) {
        dm(i) = tensor.double_data(i);
      }
    } else {
      casadi_error("Tensor has neither double_data nor raw_data");
    }

    return dm;
  }

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

      if (op_type == "Add") {
        casadi_assert(node_inputs.size() >= 2,
                      "Add operation requires 2 inputs");
        output = node_inputs[0] + node_inputs[1];

      } else if (op_type == "Mul") {
        casadi_assert(node_inputs.size() >= 2,
                      "Mul operation requires 2 inputs");
        output = node_inputs[0] * node_inputs[1];

      } else if (op_type == "Sin") {
        casadi_assert(node_inputs.size() >= 1,
                      "Sin operation requires 1 input");
        output = sin(node_inputs[0]);

      } else if (op_type == "Identity") {
        casadi_assert(node_inputs.size() >= 1,
                      "Identity operation requires 1 input");
        output = node_inputs[0];

      } else if (op_type == "Sub") {
        casadi_assert(node_inputs.size() >= 2,
                      "Sub operation requires 2 inputs");
        output = node_inputs[0] - node_inputs[1];

      } else if (op_type == "Div") {
        casadi_assert(node_inputs.size() >= 2,
                      "Div operation requires 2 inputs");
        output = node_inputs[0] / node_inputs[1];

      } else if (op_type == "Cos") {
        casadi_assert(node_inputs.size() >= 1,
                      "Cos operation requires 1 input");
        output = cos(node_inputs[0]);

      } else if (op_type == "Tan") {
        casadi_assert(node_inputs.size() >= 1,
                      "Tan operation requires 1 input");
        output = tan(node_inputs[0]);

      } else if (op_type == "Exp") {
        casadi_assert(node_inputs.size() >= 1,
                      "Exp operation requires 1 input");
        output = exp(node_inputs[0]);

      } else if (op_type == "Log") {
        casadi_assert(node_inputs.size() >= 1,
                      "Log operation requires 1 input");
        output = log(node_inputs[0]);

      } else if (op_type == "Sqrt") {
        casadi_assert(node_inputs.size() >= 1,
                      "Sqrt operation requires 1 input");
        output = sqrt(node_inputs[0]);

      } else if (op_type == "Neg") {
        casadi_assert(node_inputs.size() >= 1,
                      "Neg operation requires 1 input");
        output = -node_inputs[0];

      } else if (op_type == "Tanh") {
        casadi_assert(node_inputs.size() >= 1,
                      "Tanh operation requires 1 input");
        output = tanh(node_inputs[0]);

      } else if (op_type == "Constant") {
        // Extract constant from node attributes
        if (verbose_) {
          uout() << "    Constant node has " << node.attribute_size()
                 << " attributes" << std::endl;
        }

        casadi_assert(node.attribute_size() > 0,
                      "Constant node must have attributes");

        // Find the 'value' attribute
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
        // Unsupported operation - warn and skip
        casadi_warning("ONNX import: unsupported operation '" + op_type +
                       "' at node " + std::to_string(i) + ", skipping");
        continue;  // Skip this node
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

  void OnnxTranslator::save(const std::string& filename) {
    casadi_assert(has_model_, "No ONNX model loaded. Call load() first.");

    // Write ONNX file
    std::ofstream output(filename, std::ios::binary);
    casadi_assert(output.good(), "Failed to open output file: " + filename);

    casadi_assert(model_.SerializeToOstream(&output),
                  "Failed to serialize ONNX model to file: " + filename);

    output.close();

    if (verbose_) {
      uout() << "Saved ONNX model to: " << filename << std::endl;
    }
  }

} // namespace casadi
