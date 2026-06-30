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


#include "graph_model_internal.hpp"

namespace casadi {

  // ---------- public GraphModel handle ----------

  GraphModel::GraphModel() {
  }

  GraphModel::GraphModel(const std::string& format, const std::vector<uint8_t>& model_data,
                         const Dict& opts) {
    own(GraphModelInternal::getPlugin(format).creator(model_data));
    (*this)->construct(opts);
  }

  GraphModelInternal* GraphModel::operator->() {
    return static_cast<GraphModelInternal*>(SharedObject::operator->());
  }
  const GraphModelInternal* GraphModel::operator->() const {
    return static_cast<const GraphModelInternal*>(SharedObject::operator->());
  }
  GraphModelInternal* GraphModel::get() const {
    return static_cast<GraphModelInternal*>(SharedObject::get());
  }

  std::string GraphModel::format() const { return (*this)->plugin_name(); }

  void GraphModel::fill_metadata(GraphBuilderInternal& gb) const {
    (*this)->fill_metadata(gb);
  }
  Function GraphModel::import_symbolic(const GraphBuilderInternal& gb,
      const std::string& name) const {
    // Shallow const: the backend engine mutates internally, but importing does not
    // change the model handle's logical state.
    return get()->import_symbolic(gb, name);
  }
  std::vector<uint8_t> GraphModel::export_symbolic(const Function& f, const Dict& opts) {
    return (*this)->export_symbolic(f, opts);
  }
  const std::vector<uint8_t>& GraphModel::model_data() const { return (*this)->model_data(); }

  bool GraphModel::has_plugin(const std::string& name) {
    return GraphModelInternal::has_plugin(name);
  }
  void GraphModel::load_plugin(const std::string& name) {
    GraphModelInternal::load_plugin(name);
  }
  std::string GraphModel::doc(const std::string& name) {
    return GraphModelInternal::getPlugin(name).doc;
  }

  bool has_graphmodel(const std::string& name) { return GraphModel::has_plugin(name); }
  void load_graphmodel(const std::string& name) { GraphModel::load_plugin(name); }
  std::string doc_graphmodel(const std::string& name) { return GraphModel::doc(name); }

  // ---------- GraphModelInternal ----------

  GraphModelInternal::GraphModelInternal(const std::vector<uint8_t>& model_data)
    : model_data_(model_data), verbose_(false) {
  }

  GraphModelInternal::~GraphModelInternal() {
  }

  const Options GraphModelInternal::options_
  = {{},
     {{"verbose",
       {OT_BOOL, "Verbose evaluation -- for debugging"}}}
  };

  void GraphModelInternal::construct(const Dict& opts) {
    for (auto&& op : opts) {
      if (op.first == "verbose") verbose_ = op.second;
    }
    init(opts);
  }

  void GraphModelInternal::init(const Dict& opts) {
  }

  void GraphModelInternal::disp(std::ostream& stream, bool more) const {
    stream << "GraphModel(" << plugin_name() << ")";
  }

  std::map<std::string, GraphModelInternal::Plugin> GraphModelInternal::solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::mutex GraphModelInternal::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

  const std::string GraphModelInternal::infix_ = "graphmodel";

} // namespace casadi
