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


#ifndef CASADI_ONNX_TRANSLATOR_HPP
#define CASADI_ONNX_TRANSLATOR_HPP

#include <casadi/core/translator_internal.hpp>
#include <casadi/interfaces/onnx/casadi_translator_onnx_export.h>

#define ONNX_ML 1
#define ONNX_NAMESPACE onnx
#include <onnx/onnx_pb.h>

#include <map>
#include <string>

/// \cond INTERNAL
namespace casadi {

  /** \brief ONNX Translator

      Imports and exports ONNX computational graphs

      \author Joris Gillis
      \date 2025

      \identifier{onnx_translator} */
  class CASADI_EXPORT OnnxTranslator : public TranslatorInternal {
  public:
    /// Constructor
    OnnxTranslator();

    /// Destructor
    ~OnnxTranslator() override;

    /** \brief Get type name

        \identifier{onnx_translator_type} */
    std::string class_name() const override { return "OnnxTranslator";}

    /// Query plugin name
    const char* plugin_name() const override { return "onnx";}

    ///@{
    /** \brief Options

        \identifier{onnx_translator_options} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize

        \identifier{onnx_translator_init} */
    void init(const Dict& opts) override;

    /** \brief Load a graph from ONNX file

        \identifier{onnx_translator_load_file} */
    void load(const std::string& filename) override;

    /** \brief Load a CasADi Function and convert to ONNX representation

        \identifier{onnx_translator_load_function} */
    void load(const Function& f) override;

    /** \brief Set dimension for a symbolic variable

        \identifier{onnx_translator_set_dimension} */
    void set_dimension(const std::string& name, casadi_int dim) override;

    /** \brief Create a CasADi Function from the loaded ONNX graph

        \identifier{onnx_translator_create} */
    Function create(const std::string& name) override;

    /** \brief Save the loaded graph to ONNX file

        \identifier{onnx_translator_save} */
    void save(const std::string& filename) override;

    /// Creator function for plugin
    static TranslatorInternal* creator() {
      return new OnnxTranslator();
    }

    /// A documentation string
    static const std::string meta_doc;

  protected:
    /// ONNX model protocol buffer
    onnx::ModelProto model_;

    /// Dimension overrides for symbolic dimensions
    std::map<std::string, casadi_int> dimension_overrides_;

    /// Whether a model has been loaded
    bool has_model_;

  private:
    /** \brief Get dimension value from ONNX shape

        Handles both concrete dimensions and symbolic dimensions
        using dimension_overrides_ map.

        \identifier{onnx_translator_get_dimension} */
    casadi_int get_dimension(const onnx::TensorShapeProto& shape, int idx) const;

    /** \brief Convert ONNX TensorProto to CasADi DM

        Extracts shape and data from ONNX tensor.
        Currently only supports DOUBLE data type.

        \identifier{onnx_translator_tensor_to_dm} */
    DM tensor_to_dm(const onnx::TensorProto& tensor) const;

    /** \brief Analyze outer scope dependencies in a subgraph

        Walks through all nodes in a subgraph and identifies
        input references that are not defined locally (i.e., must
        come from outer scope).

        \param subgraph The subgraph to analyze
        \param available_in_scope Set of variable names available in outer scope
        \return Vector of variable names referenced from outer scope

        \identifier{onnx_translator_analyze_deps} */
    std::vector<std::string> analyze_outer_scope_dependencies(
        const onnx::GraphProto& subgraph,
        const std::set<std::string>& available_in_scope) const;

    /** \brief Translate ONNX subgraph to CasADi Function

        Creates a CasADi Function from an ONNX subgraph, handling
        both explicit inputs and outer scope dependencies.

        \param subgraph The ONNX subgraph (GraphProto)
        \param function_name Name for the generated Function
        \param outer_scope_vars Map of outer scope variables
        \param outer_deps Set of outer scope dependencies to include as inputs
        \return CasADi Function representing the subgraph

        \identifier{onnx_translator_translate_subgraph} */
    Function translate_subgraph_to_function(
        const onnx::GraphProto& subgraph,
        const std::string& function_name,
        const std::map<std::string, MX>& outer_scope_vars,
        const std::vector<std::string>& outer_deps);

    /** \brief Translate ONNX Loop body to CasADi Function

        Creates a CasADi Function from an ONNX Loop body subgraph.
        Special handling for iteration_num, condition, and loop-carried deps.

        \param body The Loop body subgraph (GraphProto)
        \param function_name Name for the generated Function
        \param outer_scope_vars Map of outer scope variables
        \param outer_deps Set of outer scope dependencies to include as inputs
        \return CasADi Function representing the Loop body

        \identifier{onnx_translator_translate_loop_body} */
    Function translate_loop_body_to_function(
        const onnx::GraphProto& body,
        const std::string& function_name,
        const std::map<std::string, MX>& outer_scope_vars,
        const std::vector<std::string>& outer_deps);

    /** \brief Process a single ONNX node operation

        Executes the operation logic for a single ONNX node.
        Reusable across main graph and subgraph processing.

        \param op_type ONNX operation type (e.g., "Add", "MatMul")
        \param node The ONNX NodeProto
        \param node_inputs Vector of input MX expressions
        \return MX output expression

        \identifier{onnx_translator_process_node_operation} */
    MX process_node_operation(
        const std::string& op_type,
        const onnx::NodeProto& node,
        const std::vector<MX>& node_inputs);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_ONNX_TRANSLATOR_HPP
