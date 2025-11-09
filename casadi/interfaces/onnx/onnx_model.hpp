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


#ifndef CASADI_ONNX_MODEL_HPP
#define CASADI_ONNX_MODEL_HPP

#include <casadi/core/graph_model_internal.hpp>
#include <casadi/core/mx.hpp>
#include <casadi/interfaces/onnx/casadi_graphmodel_onnx_export.h>

#define ONNX_ML 1
#define ONNX_NAMESPACE onnx
#include <onnx/onnx_pb.h>

#include <functional>
#include <map>
#include <set>
#include <string>

/// \cond INTERNAL
namespace casadi {

  class GraphBuilderInternal;

  /** \brief Callback type for adding nodes to a container (GraphProto or FunctionProto) */
  using AddNodeFn = std::function<onnx::NodeProto*()>;

  /** \brief \pluginbrief{GraphModel,onnx}

      The ONNX backend for GraphModel and the CasADi <-> ONNX graph engine in one class:
      parses the protobuf model, fills a GraphBuilder's tensor metadata, and performs
      symbolic import (graph -> MXFunction) and export (Function -> ONNX bytes). The only
      class that links onnx/protobuf.

      \author Joris Gillis
      \date 2026

      */
  class CASADI_GRAPHMODEL_ONNX_EXPORT Onnx : public GraphModelInternal {
  public:
    explicit Onnx(const std::vector<uint8_t>& model_data);
    ~Onnx() override;

    /// Plugin factory
    static GraphModelInternal* creator(const std::vector<uint8_t>& model_data) {
      return new Onnx(model_data);
    }

    const char* plugin_name() const override { return "onnx"; }
    std::string class_name() const override { return "Onnx"; }

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_; }
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    // ---- GraphModel interface ----
    void fill_metadata(GraphBuilderInternal& gb) const override;
    Function import_symbolic(const GraphBuilderInternal& gb, const std::string& name) override;
    std::vector<uint8_t> export_symbolic(const Function& f, const Dict& opts) override;

    /// Documentation
    static const std::string meta_doc;

    // ---- symbolic graph engine ----

    /** \brief Load a CasADi Function and convert to the ONNX representation */
    void load(const Function& f);

    /** \brief Set dimension for a symbolic variable */
    void set_dimension(const std::string& name, casadi_int dim);

    /** \brief Create a CasADi Function from the loaded ONNX graph */
    Function create(const std::string& name);

    /** \brief Load a graph from serialized ONNX bytes */
    void load_bytes(const std::vector<uint8_t>& data);

    /** \brief Serialize the loaded graph/model to ONNX bytes */
    std::vector<uint8_t> save_bytes() const;

  protected:
    /// ONNX model protocol buffer
    onnx::ModelProto model_;

    /// Dimension overrides for symbolic dimensions
    std::map<std::string, casadi_int> dimension_overrides_;

    /// Whether a model has been loaded
    bool has_model_ = false;

    /// Track which functions have been exported as FunctionProto
    std::set<std::string> exported_functions_;

    /// Real type for exported tensors: "double" (default) or "float"
    std::string casadi_real_ = "double";

  private:
    // Import helpers (graph -> MX), operating on this translator
    void process_graph_initializers(
        const onnx::GraphProto& graph,
        std::map<std::string, MX>& value_map,
        bool verbose) const;

    void process_graph_inputs(
        const onnx::GraphProto& graph,
        std::map<std::string, MX>& value_map,
        std::vector<MX>& func_inputs,
        std::vector<std::string>& input_names,
        bool verbose) const;

    void process_graph_nodes(
        const onnx::GraphProto& graph,
        std::map<std::string, MX>& value_map,
        bool verbose);

    void collect_graph_outputs(
        const onnx::GraphProto& graph,
        const std::map<std::string, MX>& value_map,
        std::vector<MX>& func_outputs,
        std::vector<std::string>& output_names,
        bool verbose) const;

    // Look up a model-level FunctionProto by name+domain, or nullptr if absent
    const onnx::FunctionProto* find_function(const std::string& name,
                                             const std::string& domain) const;

    /** \brief Get dimension value from ONNX shape

        Handles both concrete dimensions and symbolic dimensions
        using dimension_overrides_ map. */
    casadi_int get_dimension(const onnx::TensorShapeProto& shape, int idx) const;

    /** \brief Convert an ONNX TensorProto to a CasADi DM (DOUBLE/FLOAT/INT32/INT64/BOOL -> double) */
    DM tensor_to_dm(const onnx::TensorProto& tensor) const;

    /** \brief Convert an ONNX SparseTensorProto (COO) to a sparse CasADi DM

        Reads values [NNZ] + indices ([NNZ,2] coordinates or [NNZ] row-major linearized) and
        rebuilds the DM with DM::triplet (order-agnostic, so no reordering is needed inbound). */
    DM sparse_tensor_to_dm(const onnx::SparseTensorProto& st) const;

    /** \brief Convert a single ONNX node to an MX expression (shared by main graph and subgraphs)

        \param op_type ONNX operation type (e.g., "Add", "MatMul")
        \param node The ONNX NodeProto
        \param node_inputs Input MX expressions
        \return MX output expression */
    MX process_node_operation(
        const std::string& op_type,
        const onnx::NodeProto& node,
        const std::vector<MX>& node_inputs);

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

    /** \brief Check if Function is a reduce-map (map with reduce_in/reduce_out)

        A reduce-map flattens to an MXFunction wrapping a single MapSum, so it is
        detected by the presence of a MapSum sub-function rather than by class name. */
    bool is_reduce_map_function(const Function& f) const;

    /** \brief Error out if a called function is an unsupported control-flow construct */
    void assert_not_control_flow(const Function& called_func) const;

    /** \brief Emit a node calling a CasADi Function (exporting it as a local function once).
     *
     * Container is onnx::GraphProto or onnx::FunctionProto; outputs are named
     * out_prefix + index and recorded in work_to_onnx. */
    template<typename Container>
    void export_call(Container* container, const Function& called_func,
                     const std::vector<casadi_int>& i_vec,
                     const std::vector<casadi_int>& o,
                     std::map<casadi_int, std::string>& work_to_onnx,
                     const std::string& out_prefix);

    /** \brief Emit a Scan node for a CasADi Map (independent iteration over columns) */
    template<typename Container>
    void export_map(Container* container, const Function& map_fn,
                    const std::vector<casadi_int>& i_vec,
                    const std::vector<casadi_int>& o,
                    std::map<casadi_int, std::string>& work_to_onnx,
                    const std::string& out_prefix);

    /** \brief Build the Scan body subgraph for a Map's base function.
     *
     * Scan slices the iteration axis off the 3-D scan inputs, so the body runs the base on
     * a clean 2-D (rows, c) block per iteration. */
    onnx::GraphProto build_scan_body(const Function& base);

    /** \brief Emit a Scan node for a CasADi reduce-map.
     *
     * reduce_in inputs become outer-scope captures (broadcast every iteration); reduce_out
     * outputs become Scan state-variable accumulators; the rest are scanned/concatenated. */
    template<typename Container>
    void export_reduce_map(Container* container, const Function& wrapper,
                           const std::vector<casadi_int>& i_vec,
                           const std::vector<casadi_int>& o,
                           std::map<casadi_int, std::string>& work_to_onnx,
                           const std::string& out_prefix);

    /** \brief Build the Scan body for a reduce-map: state accumulators + captured inputs. */
    onnx::GraphProto build_reduce_scan_body(const Function& base,
                                            const std::vector<bool>& reduce_in,
                                            const std::vector<bool>& reduce_out,
                                            const std::vector<std::string>& capture_names);

    /** \brief Reconstruct a CasADi Function from an exported FunctionProto, given input shapes */
    Function function_from_function_proto(
        const onnx::FunctionProto& fp,
        const std::vector<std::pair<casadi_int, casadi_int>>& in_shapes,
        const std::string& name);

    /** \brief Emit a Constant(shape) + Reshape(data, shape) into a graph or function */
    template<typename Container>
    void emit_reshape(Container* container, const std::string& data,
                      const std::vector<casadi_int>& shape,
                      const std::string& output, const std::string& shape_name);

    /** \brief Emit a CasADi-column-major reshape (Transpose . Reshape(reverse) . Transpose) */
    template<typename Container>
    void colmajor_reshape_into(Container* container, const std::string& data,
                               const std::vector<casadi_int>& dims,
                               const std::string& output, const std::string& uniq);

    /** \brief Connect a result to a declared output: a folded reshape (source shape differs
        from the output shape) emits a colmajor reshape, otherwise a plain Identity */
    template<typename Container>
    void emit_output_node(Container* container, const std::string& data,
                          casadi_int src_rows, casadi_int src_cols, const Sparsity& out_sp,
                          const std::string& output, const std::string& uniq);

    /** \brief Build a CasADi Function from an ONNX graph (the whole model or a subgraph) */
    Function function_from_graph(const onnx::GraphProto& graph, const std::string& name);

    /** \brief Emit an If node for a CasADi Switch (2-way if_else only) */
    template<typename Container>
    void export_if(Container* container, const Function& switch_fn,
                   const std::vector<casadi_int>& i_vec,
                   const std::vector<casadi_int>& o,
                   std::map<casadi_int, std::string>& work_to_onnx,
                   const std::string& out_prefix);

    /** \brief Build an If branch subgraph: no formal inputs, captures arg_names by reference */
    onnx::GraphProto build_if_branch(const Function& f,
                                     const std::vector<std::string>& arg_names,
                                     const std::string& prefix);

    /** \brief Evaluate a captured subgraph (an If branch) against an outer name scope */
    std::vector<MX> eval_captured_subgraph(const onnx::GraphProto& graph,
                                           std::map<std::string, MX> scope);

    // --- Export helpers that depend on configuration (the real type) ---

    /** \brief ONNX element type for real-valued tensors, per the casadi_real option */
    onnx::TensorProto::DataType real_type() const {
      return casadi_real_ == "float" ? onnx::TensorProto::FLOAT : onnx::TensorProto::DOUBLE;
    }

    /** \brief Set and validate the casadi_real option ("double" or "float") */
    void set_casadi_real(const std::string& v) {
      casadi_assert(v == "double" || v == "float",
                    "casadi_real must be \"double\" or \"float\", got \"" + v + "\".");
      casadi_real_ = v;
    }

    /** \brief Declare a value as a real tensor of the given sparsity's shape */
    void set_real_tensor_type(onnx::ValueInfoProto* value, const Sparsity& sp);

    /** \brief Add graph inputs/outputs to an ONNX graph */
    void add_graph_inputs(onnx::GraphProto* graph, const Function& f,
                          const std::string& name_prefix = "");
    void add_graph_outputs(onnx::GraphProto* graph, const Function& f,
                           const std::string& name_prefix = "");

    /** \brief Add a Constant node holding a real tensor (empty dims => scalar) */
    void add_real_constant(AddNodeFn add_node, const std::string& name,
                           const std::vector<double>& data,
                           const std::vector<casadi_int>& dims = {});

    /** \brief Assemble a block matrix (diagcat / horzcat output) the pseudo-dense way: `Pad` each
        DENSE block to the output shape at its `(row_off,col_off)` offset, then `Sum`. No nonzero
        extraction -- a triangular block is just a dense block with a structural zero. The assembly
        imports as a dense block; callers restore a non-dense pattern via emit_sparsity_restore.
        Used by OP_DIAGCAT and the multi-segment output. */
    void emit_blockdiag(AddNodeFn add_node, const std::vector<std::string>& names,
                        const std::vector<casadi_int>& row_off,
                        const std::vector<casadi_int>& col_off,
                        const std::vector<casadi_int>& br, const std::vector<casadi_int>& bc,
                        casadi_int R, casadi_int C, const std::string& output,
                        const std::string& uniq);

    /** \brief Emit the general "rearrange nonzeros between two sparsities" envelope.

        Built on the poorest form shared by the whole GetNonzeros/Project family: `idx[k]` is the
        input-nonzero index feeding output-nonzero k, or -1 = fill-with-zero. Flattens `data` to a
        row, gathers (zero-filling -1 via an appended 0), and -- when `out_sp` is sparse -- scatters
        into its dense frame; then reshapes to the (transpose-rep) output. Used by OP_GETNONZEROS's
        fallback and OP_PROJECT/densify. `idx` is taken by value (remapped in place). */
    void emit_nonzero_remap(AddNodeFn add_node, const std::string& data,
                            const Sparsity& sp_in, const Sparsity& out_sp,
                            std::vector<casadi_int> idx, const std::string& uniq,
                            const std::string& node_output);

    /** \brief Plant a sparsity "seed" so an op's dense ONNX value imports to exactly target_sp.

        ONNX value-flow is dense; this emits the minimal node(s) (an Add/Mul against a
        sparse_value or dense-zero Constant) so that on import CasADi natively re-propagates
        `target_sp` from the seed. Recipes (`value_sp` is the value's own current import pattern):
         - value_sp == target_sp: no node, returns `value` unchanged (passthrough).
         - target_sp dense: Add(value, dense_zero) -> dense.
         - value_sp subset of target_sp (only drops entries): Add(value, zeros_sparse(target_sp)).
         - value_sp superset of target_sp (only adds entries, e.g. narrowing project):
           Mul(value, ones_sparse(target_sp)).
         - general: Mul(value, ones_sparse(target_sp)) then Add(_, zeros_sparse(target_sp)).
        The FINAL emitted node writes directly to `final_output` (no trailing Identity); the
        passthrough case returns `value`. Returns the name of the seeded value. */
    std::string emit_sparsity_restore(AddNodeFn add_node, const std::string& value,
                                      const Sparsity& value_sp, const Sparsity& target_sp,
                                      const std::string& uniq, const std::string& final_output);

    /** \brief Add a Constant node holding a SPARSE real tensor via the sparse_value attribute

        Stores the nonzeros as an ONNX SparseTensorProto (COO): values [NNZ], indices [NNZ,2]
        as explicit (row, col) coordinates in ONNX's required ascending row-major order, and
        the dense shape. CasADi stores nonzeros column-major; emitting them in the transpose's
        natural (column-major) order yields the original's row-major order, so no sort is needed. */
    void add_sparse_constant(AddNodeFn add_node, const std::string& name, const DM& dm);

    /** \brief Fill a SparseTensorProto with a DM's transpose-rep COO (named `name`). Shared by the
        Constant-node sparse_value path and the graph-level sparse_initializer (input overlay). */
    void fill_sparse_tensor(onnx::SparseTensorProto* st, const std::string& name,
                            const DM& dm) const;

    /** \brief Input sparsity overlay: a non-dense input's pattern is carried as a graph
        sparse_initializer (all-zero values) sharing the input's NAME -- a standard-ONNX optional
        input with a sparse zero default. Returns the recovered CasADi Sparsity, or null if absent. */
    Sparsity input_pattern(const onnx::GraphProto& graph, const std::string& name) const;

    /** \brief Process a CasADi operation into ONNX node(s); returns false if unhandled */
    bool process_operation(AddNodeFn add_node, const Function& f, casadi_int op, casadi_int k,
                           const std::vector<casadi_int>& i_vec,
                           const std::vector<casadi_int>& o_vec,
                           std::map<casadi_int, std::string>& work_to_onnx,
                           const std::string& node_output);
    bool process_operation(onnx::GraphProto* graph, const Function& f, casadi_int op, casadi_int k,
                           const std::vector<casadi_int>& i_vec,
                           const std::vector<casadi_int>& o_vec,
                           std::map<casadi_int, std::string>& work_to_onnx,
                           const std::string& node_output);
  };

  // ========== Export Helper Functions ==========

  /** \brief Function input/output name, or a generated fallback when unnamed */
  std::string onnx_input_name(const Function& f, casadi_int i);
  std::string onnx_output_name(const Function& f, casadi_int i);

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

  /** \brief Read an integer node attribute by name, or default_value if absent */
  casadi_int get_int_attribute(const onnx::NodeProto& node, const std::string& name,
                               casadi_int default_value);

  /** \brief Read a float node attribute by name, or default_value if absent */
  double get_float_attribute(const onnx::NodeProto& node, const std::string& name,
                             double default_value);

  /** \brief Read a string node attribute by name, or "" if absent */
  std::string get_string_attribute(const onnx::NodeProto& node, const std::string& name);

  /** \brief Read a subgraph (GraphProto) node attribute by name, or nullptr if absent */
  const onnx::GraphProto* get_graph_attribute(const onnx::NodeProto& node,
                                              const std::string& name);

  /** \brief Add an integer attribute (e.g. axis) to a node */
  void add_int_attribute(onnx::NodeProto* node, const std::string& name, casadi_int value);

  /** \brief Add an integer-list attribute (e.g. scan_input_axes) to a node */
  void add_ints_attribute(onnx::NodeProto* node, const std::string& name,
                          const std::vector<casadi_int>& values);

  /** \brief Add a Constant node holding an INT64 tensor to a graph (default shape: 1-D) */
  void add_int_constant(onnx::GraphProto* graph, const std::string& name,
                        const std::vector<casadi_int>& data, std::vector<casadi_int> dims = {});

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

} // namespace casadi
/// \endcond

#endif // CASADI_ONNX_MODEL_HPP
