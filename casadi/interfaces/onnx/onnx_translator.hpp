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

#include <functional>
#include <map>
#include <set>
#include <string>

/// \cond INTERNAL
namespace casadi {

  /** \brief ONNX Translator

      Imports and exports ONNX computational graphs

      \author Joris Gillis
      \date 2025

  */
  class CASADI_TRANSLATOR_ONNX_EXPORT OnnxTranslator : public TranslatorInternal {
  public:
    /// Constructor
    OnnxTranslator();

    /// Destructor
    ~OnnxTranslator() override;

    /** \brief Get type name */
    std::string class_name() const override { return "OnnxTranslator";}

    /// Query plugin name
    const char* plugin_name() const override { return "onnx";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    /** \brief Load a graph from ONNX file */
    void load(const std::string& filename) override;

    /** \brief Load a CasADi Function and convert to ONNX representation */
    void load(const Function& f) override;

    /** \brief Set dimension for a symbolic variable */
    void set_dimension(const std::string& name, casadi_int dim) override;

    /** \brief Create a CasADi Function from the loaded ONNX graph */
    Function create(const std::string& name) override;

    /** \brief Save the loaded graph to ONNX file */
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

    /// Track which functions have been exported as FunctionProto
    std::set<std::string> exported_functions_;

  private:
    // Friend declarations for import helper functions
    friend void process_graph_initializers(
        const onnx::GraphProto& graph,
        std::map<std::string, MX>& value_map,
        const OnnxTranslator& translator,
        bool verbose);

    friend void process_graph_inputs(
        const onnx::GraphProto& graph,
        std::map<std::string, MX>& value_map,
        std::vector<MX>& func_inputs,
        std::vector<std::string>& input_names,
        const OnnxTranslator& translator,
        bool verbose);

    friend void process_graph_nodes(
        const onnx::GraphProto& graph,
        std::map<std::string, MX>& value_map,
        OnnxTranslator& translator,
        bool verbose,
        bool allow_control_flow);

    friend void collect_graph_outputs(
        const onnx::GraphProto& graph,
        const std::map<std::string, MX>& value_map,
        std::vector<MX>& func_outputs,
        std::vector<std::string>& output_names,
        bool verbose);

    /** \brief Get dimension value from ONNX shape

        Handles both concrete dimensions and symbolic dimensions
        using dimension_overrides_ map. */
    casadi_int get_dimension(const onnx::TensorShapeProto& shape, int idx) const;

    /** \brief Convert ONNX TensorProto to CasADi DM

        Extracts shape and data from ONNX tensor.
        Currently only supports DOUBLE data type. */
    DM tensor_to_dm(const onnx::TensorProto& tensor) const;

    /** \brief Analyze outer scope dependencies in a subgraph

        Walks through all nodes in a subgraph and identifies
        input references that are not defined locally (i.e., must
        come from outer scope).

        \param subgraph The subgraph to analyze
        \param available_in_scope Set of variable names available in outer scope
        \return Vector of variable names referenced from outer scope

        */
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

        */
    Function translate_subgraph_to_function(
        const onnx::GraphProto& subgraph,
        const std::string& function_name,
        const std::map<std::string, MX>& outer_scope_vars,
        const std::vector<std::string>& outer_deps);

    /** \brief Process a single ONNX node operation

        Executes the operation logic for a single ONNX node.
        Reusable across main graph and subgraph processing.

        \param op_type ONNX operation type (e.g., "Add", "MatMul")
        \param node The ONNX NodeProto
        \param node_inputs Vector of input MX expressions
        \return MX output expression */
    MX process_node_operation(
        const std::string& op_type,
        const onnx::NodeProto& node,
        const std::vector<MX>& node_inputs);

    /** \brief Convert a CasADi Function to ONNX GraphProto (for subgraphs)

        Recursively converts a CasADi Function into an ONNX GraphProto.
        Used for exporting control flow operators (If, Loop, Scan).

        \param f The CasADi Function to convert
        \param graph_name Name for the ONNX graph
        \param outer_scope_inputs Optional vector of outer scope variable names to use
               instead of adding graph inputs. Used for If branch subgraphs that
               reference variables from the enclosing scope.
        \return Pointer to created GraphProto */
    onnx::GraphProto* function_to_graph(
        const Function& f,
        const std::string& graph_name,
        const std::vector<std::string>& outer_scope_inputs = {});

    /** \brief Convert a CasADi Function to ONNX FunctionProto

        Converts a CasADi Function into an ONNX FunctionProto for use as
        a local function that can be called from the main graph or other
        functions. Recursively handles nested function calls.

        \param f The CasADi Function to convert
        \param domain The ONNX domain for this function (e.g., "casadi")
        \return Pointer to created FunctionProto */
    onnx::FunctionProto* function_to_function_proto(
        const Function& f,
        const std::string& domain);

    /** \brief Check if Function is an if_else function

        \param f The Function to check
        \return True if this is an if_else function */
    bool is_if_else_function(const Function& f) const;

    /** \brief Check if Function is a mapaccum function

        \param f The Function to check
        \return True if this is a mapaccum function */
    bool is_mapaccum_function(const Function& f) const;

    /** \brief Check if Function is a map function

        \param f The Function to check
        \return True if this is a map function

    */
    bool is_map_function(const Function& f) const;

    /** \brief Create ONNX node from CasADi operation

        Shared helper to convert a CasADi operation to an ONNX node.
        Used by both main export and subgraph export.

        \param graph The ONNX graph to add nodes to
        \param op CasADi operation code
        \param node_output Output name for the node
        \param i_vec Input work vector indices
        \param work_to_onnx Map from work indices to ONNX names
        \param f The Function being exported
        \param k Instruction index
        \return True if operation was handled, false otherwise */
    bool create_onnx_node_from_op(
        onnx::GraphProto* graph,
        casadi_int op,
        const std::string& node_output,
        const std::vector<casadi_int>& i_vec,
        const std::vector<casadi_int>& o_vec,
        std::map<casadi_int, std::string>& work_to_onnx,
        const Function& f,
        casadi_int k);
  };

  // ========== Import Helper Functions (onnx_import.cpp) ==========

  /** \brief Process ONNX initializers (constants) */
  void process_graph_initializers(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      const OnnxTranslator& translator,
      bool verbose);

  /** \brief Create MX symbols for graph inputs */
  void process_graph_inputs(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      std::vector<MX>& func_inputs,
      std::vector<std::string>& input_names,
      const OnnxTranslator& translator,
      bool verbose);

  /** \brief Process all nodes in the graph */
  void process_graph_nodes(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      OnnxTranslator& translator,
      bool verbose,
      bool allow_control_flow = true);

  /** \brief Collect graph outputs from value_map */
  void collect_graph_outputs(
      const onnx::GraphProto& graph,
      const std::map<std::string, MX>& value_map,
      std::vector<MX>& func_outputs,
      std::vector<std::string>& output_names,
      bool verbose);

  // ========== Export Helper Functions ==========

  /** \brief Add graph inputs to ONNX graph (defined in onnx_export.cpp) */
  void add_graph_inputs(onnx::GraphProto* graph, const Function& f,
                        const std::string& name_prefix = "");

  /** \brief Add graph outputs to ONNX graph (defined in onnx_export.cpp) */
  void add_graph_outputs(onnx::GraphProto* graph, const Function& f,
                         const std::string& name_prefix = "");

  /** \brief Operation mapping between CasADi and ONNX
   *
   * Single source of truth for all simple operations that have
   * direct CasADi ↔ ONNX equivalents. */
  struct OpMapping {
    casadi_int casadi_op;  ///< CasADi operation code
    const char* onnx_name; ///< ONNX operation name
    int arity;             ///< Number of inputs (1=unary, 2=binary)
  };

  /** \brief Lookup operation mapping by CasADi opcode (for export)
   *
   * Returns nullptr if not a simple mapped operation. */
  const OpMapping* get_op_mapping(casadi_int op);

  /** \brief Lookup operation mapping by ONNX name (for import)
   *
   * Returns nullptr if not a simple mapped operation. */
  const OpMapping* get_op_mapping_by_name(const std::string& onnx_name);

  /** \brief Callback type for adding nodes to a container (GraphProto or FunctionProto) */
  using AddNodeFn = std::function<onnx::NodeProto*()>;

  /** \brief Create binary operation ONNX node (callback-based) */
  onnx::NodeProto* create_binary_node(
      AddNodeFn add_node,
      const std::string& op_type,
      const std::string& input1,
      const std::string& input2,
      const std::string& output);

  /** \brief Create unary operation ONNX node (callback-based) */
  onnx::NodeProto* create_unary_node(
      AddNodeFn add_node,
      const std::string& op_type,
      const std::string& input,
      const std::string& output);

  /** \brief Process a CasADi operation and create corresponding ONNX node (callback-based) */
  bool process_operation(
      AddNodeFn add_node,
      const Function& f,
      casadi_int op,
      casadi_int k,
      const std::vector<casadi_int>& i_vec,
      const std::vector<casadi_int>& o_vec,
      std::map<casadi_int, std::string>& work_to_onnx,
      const std::string& node_output);

  /** \brief Create binary operation ONNX node (GraphProto convenience wrapper) */
  onnx::NodeProto* create_binary_node(
      onnx::GraphProto* graph,
      const std::string& op_type,
      const std::string& input1,
      const std::string& input2,
      const std::string& output);

  /** \brief Create unary operation ONNX node (GraphProto convenience wrapper) */
  onnx::NodeProto* create_unary_node(
      onnx::GraphProto* graph,
      const std::string& op_type,
      const std::string& input,
      const std::string& output);

  /** \brief Process a CasADi operation and create corresponding ONNX node (GraphProto convenience wrapper) */
  bool process_operation(
      onnx::GraphProto* graph,
      const Function& f,
      casadi_int op,
      casadi_int k,
      const std::vector<casadi_int>& i_vec,
      const std::vector<casadi_int>& o_vec,
      std::map<casadi_int, std::string>& work_to_onnx,
      const std::string& node_output);

} // namespace casadi
/// \endcond

#endif // CASADI_ONNX_TRANSLATOR_HPP
