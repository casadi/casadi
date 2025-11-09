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

#include "graph_builder_internal.hpp"
#include "onnx_function_impl.hpp"
#include <fstream>

namespace casadi {

  // ---------- Node ----------

  void Node::disp(std::ostream& stream, bool more) const {
    (void)more;
    stream << io << " " << name << ": " << dtype << "[";
    for (size_t k = 0; k < shape.size(); ++k) {
      if (k) stream << "x";
      if (shape[k] >= 0) stream << shape[k];
      else
        stream << (dim_params[k].empty() ? "?" : dim_params[k]);
    }
    stream << "]";
    if (baked) stream << " (baked)";
  }

  std::string Node::get_str(bool more) const {
    std::stringstream ss;
    disp(ss, more);
    return ss.str();
  }

  // Infer the model format from a file suffix
  static std::string format_from_path(const std::string& path) {
    auto dot = path.find_last_of('.');
    std::string suffix = dot == std::string::npos ? "" : path.substr(dot + 1);
    if (suffix == "onnx") return "onnx";
    casadi_error("GraphBuilder: cannot infer format from '" + path + "'");
  }

  // ---------- public GraphBuilder ----------

  GraphBuilder::GraphBuilder() {
  }

  GraphBuilder::GraphBuilder(const std::string& model_path, const Dict& opts) {
    std::ifstream file(model_path, std::ios::binary | std::ios::ate);
    casadi_assert(file.is_open(), "Cannot open model file: " + model_path);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    std::vector<uint8_t> data(static_cast<size_t>(size));
    casadi_assert(file.read(reinterpret_cast<char*>(data.data()), size),
                  "Cannot read model file: " + model_path);
    own(new GraphBuilderInternal(model_path, data, format_from_path(model_path), opts));
  }

  GraphBuilder::GraphBuilder(const Function& f, const Dict& opts) {
    own(new GraphBuilderInternal(f.name(), f, opts));
  }

  GraphBuilder::GraphBuilder(const std::string& name, const std::vector<uint8_t>& model_data,
                             const std::string& format, const Dict& opts) {
    own(new GraphBuilderInternal(name, model_data, format, opts));
  }

  GraphBuilderInternal* GraphBuilder::operator->() {
    return static_cast<GraphBuilderInternal*>(SharedObject::operator->());
  }
  const GraphBuilderInternal* GraphBuilder::operator->() const {
    return static_cast<const GraphBuilderInternal*>(SharedObject::operator->());
  }
  GraphBuilderInternal* GraphBuilder::get() const {
    return static_cast<GraphBuilderInternal*>(SharedObject::get());
  }

  const std::string& GraphBuilder::name() const {
    static std::string null = "null";
    return is_null() ? null : (*this)->name_;
  }

  casadi_int GraphBuilder::n_in() const { return (*this)->n_in(); }
  casadi_int GraphBuilder::n_out() const { return (*this)->n_out(); }
  std::vector<std::string> GraphBuilder::name_in() const { return (*this)->name_in(); }
  std::vector<std::string> GraphBuilder::name_out() const { return (*this)->name_out(); }
  std::vector<casadi_int> GraphBuilder::input_shape(const std::string& name) const {
    return (*this)->input_shape(name);
  }
  std::vector<casadi_int> GraphBuilder::output_shape(const std::string& name) const {
    return (*this)->output_shape(name);
  }
  std::vector<std::string> GraphBuilder::dynamic_params() const {
    return (*this)->dynamic_params();
  }
  Node GraphBuilder::node(const std::string& name) const { return (*this)->node(name); }
  std::vector<Node> GraphBuilder::nodes() const { return (*this)->nodes(); }

  void GraphBuilder::bind_dim(const std::string& param, casadi_int value) {
    (*this)->bind_dim(param, value);
  }
  void GraphBuilder::bind_shape(const std::string& input_name,
                                const std::vector<casadi_int>& shape) {
    (*this)->bind_shape(input_name, shape);
  }
  void GraphBuilder::set(const std::string& input_name, const std::vector<double>& value) {
    (*this)->set_value(input_name, value);
  }
  void GraphBuilder::set(const std::string& input_name, double value) {
    (*this)->set_value(input_name, std::vector<double>(1, value));
  }

  Function GraphBuilder::create(const std::string& name,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts) const {
    return (*this)->create_function(name, name_in, name_out, opts);
  }
  Function GraphBuilder::create(const std::string& name, const Dict& opts) const {
    return (*this)->create_function(name, {}, {}, opts);
  }
  void GraphBuilder::export_onnx(const std::string& filename, const Dict& opts) {
    (*this)->export_onnx(filename, opts);
  }

  // ---------- GraphBuilderInternal ----------

  GraphBuilderInternal::GraphBuilderInternal(const std::string& name,
                                             const std::vector<uint8_t>& model_data,
                                             const std::string& format, const Dict& opts)
    : name_(name), format_(format), model_data_(model_data) {
    model_ = GraphModel(format, model_data, opts);
    model_.fill_metadata(*this);
  }

  GraphBuilderInternal::GraphBuilderInternal(const std::string& name, const Function& f,
                                             const Dict& opts)
    : name_(name), format_("onnx"), fun_(f) {
    populate_from_function();
  }

  GraphBuilderInternal::~GraphBuilderInternal() {
  }

  void GraphBuilderInternal::populate_from_function() {
    nodes_.clear();
    for (casadi_int i = 0; i < fun_.n_in(); ++i) {
      Node n;
      n.name = fun_.name_in(i);
      n.io = "input";
      n.shape = {fun_.size1_in(i), fun_.size2_in(i)};
      n.dim_params = {"", ""};
      nodes_.push_back(n);
    }
    for (casadi_int i = 0; i < fun_.n_out(); ++i) {
      Node n;
      n.name = fun_.name_out(i);
      n.io = "output";
      n.shape = {fun_.size1_out(i), fun_.size2_out(i)};
      n.dim_params = {"", ""};
      nodes_.push_back(n);
    }
  }

  std::vector<casadi_int> GraphBuilderInternal::resolved_shape(const Node& n) const {
    if (n.io == "input") {
      auto ov = input_shapes_.find(n.name);
      if (ov != input_shapes_.end()) return ov->second;  // explicit override
    }
    std::vector<casadi_int> shape;
    for (size_t k = 0; k < n.shape.size(); ++k) {
      casadi_int d = n.shape[k];
      if (d < 0) {  // dynamic: bound name else default 1
        auto it = dim_bindings_.find(n.dim_params[k]);
        d = (it != dim_bindings_.end()) ? it->second : 1;
      }
      shape.push_back(d);
    }
    return shape;
  }

  const Node& GraphBuilderInternal::find(const std::string& name, const std::string& io) const {
    for (const Node& n : nodes_) if (n.io == io && n.name == name) return n;
    casadi_error("Graph tensor '" + name + "' (" + io + ") not found in model '" + name_ + "'");
  }

  casadi_int GraphBuilderInternal::n_in() const {
    casadi_int c = 0;
    for (const Node& n : nodes_) if (n.io == "input") ++c;
    return c;
  }
  casadi_int GraphBuilderInternal::n_out() const {
    casadi_int c = 0;
    for (const Node& n : nodes_) if (n.io == "output") ++c;
    return c;
  }
  std::vector<std::string> GraphBuilderInternal::name_in() const {
    std::vector<std::string> r;
    for (const Node& n : nodes_) if (n.io == "input") r.push_back(n.name);
    return r;
  }
  std::vector<std::string> GraphBuilderInternal::name_out() const {
    std::vector<std::string> r;
    for (const Node& n : nodes_) if (n.io == "output") r.push_back(n.name);
    return r;
  }
  std::vector<casadi_int> GraphBuilderInternal::input_shape(const std::string& name) const {
    return find(name, "input").shape;
  }
  std::vector<casadi_int> GraphBuilderInternal::output_shape(const std::string& name) const {
    return find(name, "output").shape;
  }

  std::vector<std::string> GraphBuilderInternal::dynamic_params() const {
    std::vector<std::string> r;
    for (const Node& n : nodes_) {
      for (size_t k = 0; k < n.shape.size(); ++k) {
        if (n.shape[k] < 0 && !n.dim_params[k].empty() &&
            std::find(r.begin(), r.end(), n.dim_params[k]) == r.end()) {
          r.push_back(n.dim_params[k]);
        }
      }
    }
    return r;
  }

  Node GraphBuilderInternal::node(const std::string& name) const {
    for (const Node& n : nodes_) if (n.name == name) return n;
    casadi_error("Graph tensor '" + name + "' not found in model '" + name_ + "'");
  }

  void GraphBuilderInternal::set_value(const std::string& input_name,
                                       const std::vector<double>& value) {
    find(input_name, "input");  // validate the name
    input_values_[input_name] = value;
    for (Node& n : nodes_) if (n.io == "input" && n.name == input_name) {
      n.value = value;
      n.baked = true;
    }
  }

  void GraphBuilderInternal::bind_shape(const std::string& input_name,
                                        const std::vector<casadi_int>& shape) {
    const Node& t = find(input_name, "input");
    casadi_assert(shape.size() == t.shape.size(),
                  "bind_shape: rank mismatch for '" + input_name + "'");
    input_shapes_[input_name] = shape;
    // Pinning a named dynamic axis also binds that dim everywhere it appears (e.g. outputs)
    for (size_t k = 0; k < t.shape.size(); ++k) {
      if (t.shape[k] < 0 && !t.dim_params[k].empty()) dim_bindings_[t.dim_params[k]] = shape[k];
    }
  }

  Function GraphBuilderInternal::create_function(const std::string& name,
      const std::vector<std::string>& inputs,
      const std::vector<std::string>& outputs,
      const Dict& opts) const {
    bool symbolic = false;
    std::string backend = "ort";
    Dict o;
    for (auto&& op : opts) {
      if (op.first == "symbolic") symbolic = op.second;
      else if (op.first == "backend") backend = op.second.to_string();
      else
        o[op.first] = op.second;
    }

    if (symbolic) {
      casadi_assert(!model_.is_null(),
        "GraphBuilder: symbolic create requires a parsed model (build from a file)");
      return model_.import_symbolic(*this, name);
    }

    // Numeric path: OnnxFunction freezes a snapshot directly from this builder's config
    casadi_assert(!model_data_.empty(), "GraphBuilder: numeric create requires model bytes");
    return OnnxFunction::create(backend, name, this, inputs, outputs, o);
  }

  void GraphBuilderInternal::export_onnx(const std::string& filename, const Dict& opts) {
    std::vector<uint8_t> bytes;
    if (!fun_.is_null()) {
      GraphModel gm(format_);
      bytes = gm.export_symbolic(fun_, opts);
    } else {
      casadi_assert(!model_data_.empty(), "GraphBuilder: nothing to export");
      bytes = model_data_;
    }
    std::ofstream out(filename, std::ios::binary);
    casadi_assert(out.good(), "Cannot open output file: " + filename);
    out.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
  }

  void GraphBuilderInternal::disp(std::ostream& stream, bool more) const {
    stream << "GraphBuilder '" << name_ << "': " << n_in() << " input(s), "
           << n_out() << " output(s)";
    if (!more) return;
    stream << "\nInputs:";
    for (const Node& n : nodes_) if (n.io == "input") stream << "\n  " << n.get_str();
    stream << "\nOutputs:";
    for (const Node& n : nodes_) if (n.io == "output") stream << "\n  " << n.get_str();
    std::vector<std::string> dp = dynamic_params();
    if (!dp.empty()) {
      stream << "\nDynamic dimensions:";
      for (const std::string& p : dp) stream << " " << p;
    }
  }

} // namespace casadi
