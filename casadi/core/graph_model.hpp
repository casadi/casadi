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


#ifndef CASADI_GRAPH_MODEL_HPP
#define CASADI_GRAPH_MODEL_HPP

#include "function.hpp"

namespace casadi {

  class GraphModelInternal;
  class GraphBuilderInternal;

#ifndef SWIG
  /** \brief Format-agnostic handle to a parsed computational-graph model

      Loaded by format name as a plugin (e.g. "onnx"); the concrete backend is the only
      place a heavy dependency (protobuf) is linked. The Fmu analogue for GraphBuilder.
      Internal: owned by GraphBuilderInternal, not exposed to the scripting layer.

      \date 2026

      */
  class CASADI_EXPORT GraphModel : public SharedObject {
   public:
    /// Readable name of the class
    std::string type_name() const { return "GraphModel"; }

    /// Default constructor (null handle)
    GraphModel();

    /// Load a backend for a format, optionally with parsed model bytes
    explicit GraphModel(const std::string& format,
                        const std::vector<uint8_t>& model_data = std::vector<uint8_t>(),
                        const Dict& opts = Dict());

    /// Format name (plugin name)
    std::string format() const;

    /// Populate a builder's Node metadata from the parsed model
    void fill_metadata(GraphBuilderInternal& gb) const;

    /// Symbolic import: rebuild the graph as a CasADi Function
    Function import_symbolic(const GraphBuilderInternal& gb, const std::string& name) const;

    /// Symbolic export: serialize a CasADi Function to model bytes
    std::vector<uint8_t> export_symbolic(const Function& f, const Dict& opts = Dict());

    /// Raw model bytes
    const std::vector<uint8_t>& model_data() const;

#ifndef SWIG
    ///@{
    /// Access functions of the node
    GraphModelInternal* operator->();
    const GraphModelInternal* operator->() const;
    GraphModelInternal* get() const;
    ///@}
#endif // SWIG

    /// Check if a format plugin is available
    static bool has_plugin(const std::string& name);

    /// Explicitly load a format plugin
    static void load_plugin(const std::string& name);

    /// Get format-specific documentation
    static std::string doc(const std::string& name);
  };
#endif // SWIG

  /// Check if a graph-model format plugin is available
  CASADI_EXPORT bool has_graphmodel(const std::string& name);

  /// Explicitly load a graph-model format plugin
  CASADI_EXPORT void load_graphmodel(const std::string& name);

  /// Get format-specific documentation
  CASADI_EXPORT std::string doc_graphmodel(const std::string& name);

} // namespace casadi

#endif // CASADI_GRAPH_MODEL_HPP
