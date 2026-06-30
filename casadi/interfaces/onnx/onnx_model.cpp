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
#include <casadi/core/graph_builder_internal.hpp>

namespace casadi {

  // Defined in casadi core (onnx.cpp): human-readable name of an ONNX element-type enum
  std::string onnx_dtype_name(casadi_int elem_type);

  extern "C"
  int CASADI_GRAPHMODEL_ONNX_EXPORT
  casadi_register_graphmodel_onnx(GraphModelInternal::Plugin* plugin) {
    plugin->creator = Onnx::creator;
    plugin->name = "onnx";
    plugin->doc = Onnx::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Onnx::options_;
    return 0;
  }

  extern "C"
  void CASADI_GRAPHMODEL_ONNX_EXPORT casadi_load_graphmodel_onnx() {
    GraphModelInternal::registerPlugin(casadi_register_graphmodel_onnx);
  }

  const std::string Onnx::meta_doc =
    "ONNX backend for GraphModel: protobuf metadata, symbolic import and export.\n";

  const Options Onnx::options_
  = {{&GraphModelInternal::options_},
     {{"casadi_real",
       {OT_STRING, "Real type for exported tensors: 'double' (default) or 'float'"}}
     }
    };

  Onnx::Onnx(const std::vector<uint8_t>& model_data) : GraphModelInternal(model_data) {
    if (!model_data.empty()) load_bytes(model_data);
  }

  Onnx::~Onnx() {
  }

  void Onnx::init(const Dict& opts) {
    GraphModelInternal::init(opts);
    for (auto&& op : opts) {
      if (op.first == "casadi_real") set_casadi_real(op.second.to_string());
    }
  }

  void Onnx::load_bytes(const std::vector<uint8_t>& data) {
    casadi_assert(model_.ParseFromArray(data.data(), static_cast<int>(data.size())),
                  "Failed to parse ONNX model from memory");
    has_model_ = true;
  }

  std::vector<uint8_t> Onnx::save_bytes() const {
    casadi_assert(has_model_, "No ONNX model loaded.");
    std::string s;
    casadi_assert(model_.SerializeToString(&s), "Failed to serialize ONNX model");
    return std::vector<uint8_t>(s.begin(), s.end());
  }

  // Build a Node descriptor from a graph input/output ValueInfo
  static Node io_node(const onnx::ValueInfoProto& vi, const std::string& io) {
    Node n;
    n.name = vi.name();
    n.io = io;
    const onnx::TypeProto::Tensor& tt = vi.type().tensor_type();
    n.dtype = onnx_dtype_name(tt.elem_type());
    const onnx::TensorShapeProto& sh = tt.shape();
    for (int k = 0; k < sh.dim_size(); ++k) {
      const auto& d = sh.dim(k);
      if (d.has_dim_value()) {
        n.dimension.push_back(static_cast<casadi_int>(d.dim_value()));
        n.dim_params.push_back("");
      } else {  // symbolic or unspecified -> dynamic
        n.dimension.push_back(-1);
        n.dim_params.push_back(d.has_dim_param() ? d.dim_param() : "");
      }
    }
    return n;
  }

  void Onnx::fill_metadata(GraphBuilderInternal& gb) const {
    gb.clear_nodes();
    const onnx::GraphProto& g = model_.graph();

    // Initializers are constants, not graph inputs (mirror ORT, which excludes them)
    std::set<std::string> init_names;
    for (int i = 0; i < g.initializer_size(); ++i) init_names.insert(g.initializer(i).name());

    for (int i = 0; i < g.input_size(); ++i)
      if (!init_names.count(g.input(i).name())) gb.add_node(io_node(g.input(i), "input"));
    for (int i = 0; i < g.output_size(); ++i) gb.add_node(io_node(g.output(i), "output"));
  }

  Function Onnx::import_symbolic(const GraphBuilderInternal& gb, const std::string& name) {
    for (auto&& b : gb.dim_bindings()) set_dimension(b.first, b.second);
    return create(name);
  }

  std::vector<uint8_t> Onnx::export_symbolic(const Function& f, const Dict& opts) {
    auto it = opts.find("casadi_real");
    if (it != opts.end()) set_casadi_real(it->second.to_string());
    load(f);
    return save_bytes();
  }

} // namespace casadi
