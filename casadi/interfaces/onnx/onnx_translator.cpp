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

  Function OnnxTranslator::create(const std::string& name) {
    casadi_assert(has_model_, "No ONNX model loaded. Call load() first.");

    // TODO: Implement ONNX to CasADi Function conversion
    // This is where you'll translate the ONNX graph to CasADi MX expressions
    casadi_error("OnnxTranslator::create() not yet implemented");

    return Function();
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
