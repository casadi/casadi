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


#ifndef CASADI_GRAPH_MODEL_INTERNAL_HPP
#define CASADI_GRAPH_MODEL_INTERNAL_HPP

#include "graph_model.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  class GraphBuilderInternal;

  /** \brief Base interface for format-specific graph-model backends

      Subclasses (e.g. Onnx) own the heavy dependency and implement metadata fill,
      symbolic import and symbolic export. Registered as plugins by format name.

      \date 2026

      */
  class CASADI_EXPORT GraphModelInternal
    : public SharedObjectInternal,
      public PluginInterface<GraphModelInternal> {
   public:
    explicit GraphModelInternal(const std::vector<uint8_t>& model_data);
    ~GraphModelInternal() override;

    std::string class_name() const override { return "GraphModelInternal"; }
    void disp(std::ostream& stream, bool more) const override;

    /// Plugin creator function type
    typedef GraphModelInternal* (*Creator)(const std::vector<uint8_t>& model_data);

    /// No statically-exposed plugin functions
    struct Exposed { };

    ///@{
    /** \brief Options */
    static const Options options_;
    virtual const Options& get_options() const { return options_; }
    ///@}

    /// Prepare the backend for use
    void construct(const Dict& opts);

    /** \brief Initialize */
    virtual void init(const Dict& opts);

    /// Populate a builder's Node metadata from the parsed model
    virtual void fill_metadata(GraphBuilderInternal& gb) const = 0;

    /// Rebuild the graph as a CasADi Function (symbolic import; mutates the backend's engine)
    virtual Function import_symbolic(const GraphBuilderInternal& gb,
                                     const std::string& name) = 0;

    /// Serialize a CasADi Function as model bytes (symbolic export; mutates the backend's engine)
    virtual std::vector<uint8_t> export_symbolic(const Function& f, const Dict& opts) = 0;

    /// Raw model bytes
    const std::vector<uint8_t>& model_data() const { return model_data_; }

    /// Query plugin name
    const char* plugin_name() const override = 0;

    /// Collection of available format plugins
    static std::map<std::string, Plugin> solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /// Infix used to form plugin registration symbols (casadi_register_graphmodel_<name>)
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "graphmodel"; }

   protected:
    /// Raw model bytes (empty when constructed for export only)
    std::vector<uint8_t> model_data_;

    /// Verbose -- for debugging
    bool verbose_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_GRAPH_MODEL_INTERNAL_HPP
