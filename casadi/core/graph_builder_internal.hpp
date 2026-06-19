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


#ifndef CASADI_GRAPH_BUILDER_INTERNAL_HPP
#define CASADI_GRAPH_BUILDER_INTERNAL_HPP

#include "graph_builder.hpp"
#include "graph_model_impl.hpp"
#include "shared_object.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Metadata for one graph tensor (graph input or output)

      The format-neutral GraphBuilder analogue of DaeBuilder's Variable.

      \identifier{2ju} */
  struct CASADI_EXPORT Node
    : public SWIG_IF_ELSE(PrintableCommon, Printable<Node>) {
    std::string name;                     ///< Tensor name
    std::string io;                       ///< "input" or "output"
    std::string dtype;                    ///< Element type name (FLOAT, INT64, ...)
    std::vector<casadi_int> dimension;    ///< Declared shape, -1 for dynamic dimensions
    std::vector<std::string> dim_params;  ///< Symbolic name per axis ("" if static)
    std::vector<double> value;            ///< Baked value (optional)
    bool baked = false;                   ///< True if a fixed value was set() for this input

    std::string type_name() const { return "Node"; }
    void disp(std::ostream& stream, bool more=false) const;
    std::string get_str(bool more=false) const;
  };

  /** \brief Internal class for GraphBuilder

      Single mutable holder of graph metadata + configuration; dependency-free. Metadata is
      populated either from a parsed model (via a GraphModel backend) or from a source Function.

      \date 2026

      \identifier{2jv} */
  class CASADI_EXPORT GraphBuilderInternal : public SharedObjectInternal {
   public:
    /// Construct from parsed model bytes of a given format
    GraphBuilderInternal(const std::string& name, const std::vector<uint8_t>& model_data,
                         const std::string& format, const Dict& opts);
    /// Construct from a source Function (export lifecycle)
    GraphBuilderInternal(const std::string& name, const Function& f, const Dict& opts);
    ~GraphBuilderInternal() override;

    std::string class_name() const override { return "GraphBuilderInternal"; }
    void disp(std::ostream& stream, bool more) const override;

    casadi_int n_in() const;
    casadi_int n_out() const;
    std::vector<std::string> name_in() const;
    std::vector<std::string> name_out() const;
    std::vector<casadi_int> input_shape(const std::string& name) const;
    std::vector<casadi_int> output_shape(const std::string& name) const;
    std::vector<std::string> dynamic_params() const;
    Node node(const std::string& name) const;
    std::vector<Node> nodes() const { return nodes_; }

    void bind_dim(const std::string& param, casadi_int value) { dim_bindings_[param] = value; }
    void bind_shape(const std::string& input_name, const std::vector<casadi_int>& shape);
    void set_value(const std::string& input_name, const std::vector<double>& value);

    Function create_function(const std::string& name,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts) const;
    void export_onnx(const std::string& filename, const Dict& opts);

    /// Append a tensor descriptor (called by a backend during fill_metadata)
    void add_node(const Node& n) { nodes_.push_back(n); }
    /// Drop all tensor descriptors (called by a backend before re-filling)
    void clear_nodes() { nodes_.clear(); }

    ///@{
    /// Read-only views for backends
    const std::vector<Node>& node_list() const { return nodes_; }
    const std::map<std::string, casadi_int>& dim_bindings() const { return dim_bindings_; }
    ///@}

    /// Resolve a node's declared shape to concrete sizes: input_shapes_ override (inputs),
    /// else dynamic axes bound via dim_bindings_ (named) or defaulted to 1
    std::vector<casadi_int> resolved_shape(const Node& n) const;

    /// Locate a node by name in a given I/O role (throws if absent)
    const Node& find(const std::string& name, const std::string& io) const;

    std::string name_;
    std::string format_;
    std::vector<uint8_t> model_data_;

    /// Source Function (export lifecycle); null when built from a model
    Function fun_;
    /// Parsed model backend (import lifecycle); null when built from a Function
    GraphModel model_;

    /// Tensor metadata (inputs followed by outputs)
    std::vector<Node> nodes_;

    /// Pending configuration carried into create()
    std::map<std::string, casadi_int> dim_bindings_;
    std::map<std::string, std::vector<casadi_int>> input_shapes_;
    std::map<std::string, std::vector<double>> input_values_;

   private:
    /// Derive metadata from the source Function (names, shapes)
    void populate_from_function();
  };

} // namespace casadi

/// \endcond

#endif // CASADI_GRAPH_BUILDER_INTERNAL_HPP
