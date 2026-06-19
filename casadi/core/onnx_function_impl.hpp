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


#ifndef CASADI_ONNX_FUNCTION_IMPL_HPP
#define CASADI_ONNX_FUNCTION_IMPL_HPP

#include "onnx_function.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"
#include <set>

/// \cond INTERNAL

namespace casadi {

  class GraphBuilderInternal;

  /// Metadata for a single ONNX tensor (input or output); shape is resolved by the builder
  struct OnnxTensorInfo {
    std::string name;                    ///< ONNX tensor name
    std::vector<casadi_int> shape;       ///< Resolved shape (dynamic dims bound or set to 1)
    casadi_int elem_type;                ///< ONNX element type enum (1=float, 11=double, 7=int64)
    casadi_int numel;                    ///< Number of elements in the resolved shape
  };

  struct CASADI_EXPORT OnnxMemory : public FunctionMemory {
  };

  /// Human-readable name of an ONNX element-type enum (1=FLOAT, 11=DOUBLE, 7=INT64, ...)
  CASADI_EXPORT std::string onnx_dtype_name(casadi_int elem_type);

  /// ONNX element-type enum for a human-readable name (inverse of onnx_dtype_name; 0 if unknown)
  CASADI_EXPORT casadi_int onnx_dtype_enum(const std::string& name);

  /** \brief Black-box ONNX function base; backends (e.g. onnxruntime) derive from this

      \date 2024

      */
  class CASADI_EXPORT OnnxFunction : public FunctionInternal, public PluginInterface<OnnxFunction> {
   public:
    /// Construct by freezing a snapshot of a builder's metadata + config (exposed selection)
    OnnxFunction(const std::string& name,
                 const GraphBuilderInternal* gb,
                 const std::vector<std::string>& inputs,
                 const std::vector<std::string>& outputs);
    ~OnnxFunction() override;

    /// Plugin creator function type
    typedef OnnxFunction* (*Creator)(const std::string& name,
                                     const GraphBuilderInternal* gb,
                                     const std::vector<std::string>& inputs,
                                     const std::vector<std::string>& outputs,
                                     const Dict& opts);

    /// No statically-exposed plugin functions
    struct Exposed { };

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // The plugin fills tensor metadata before init(); the base exposes I/O from it
    size_t get_n_in() override { return in_.size(); }
    size_t get_n_out() override { return out_.size(); }
    std::string get_name_in(casadi_int i) override { return in_.at(i).name; }
    std::string get_name_out(casadi_int i) override { return out_.at(i).name; }
    Sparsity get_sparsity_in(casadi_int i) override { return tensor_sparsity(in_.at(i).shape); }
    Sparsity get_sparsity_out(casadi_int i) override { return tensor_sparsity(out_.at(i).shape); }

    /** \brief Initialize */
    void init(const Dict& opts) override;

    void* alloc_mem() const override { return new OnnxMemory(); }
    void free_mem(void* mem) const override { delete static_cast<OnnxMemory*>(mem); }

    int eval(const double** arg, double** res,
             casadi_int* iw, double* w, void* mem) const override = 0;

    bool uses_output() const override { return false; }

    // A model authored with CasADi's derivative naming can serve its own forward sensitivities:
    // get_forward re-create()s from the same model selecting the fwd_* tensors. Otherwise CasADi
    // falls back to finite differences.
    bool has_forward(casadi_int nfwd) const override;
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    bool has_reverse(casadi_int nadj) const override;
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    bool has_jacobian() const override;
    Function get_jacobian(const std::string& name,
                          const std::vector<std::string>& inames,
                          const std::vector<std::string>& onames,
                          const Dict& opts) const override;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;
    /** \brief Serialize type information */
    void serialize_type(SerializingStream &s) const override;
    /** \brief String used to identify the immediate FunctionInternal subclass */
    std::string serialize_base_function() const override { return "Onnx"; }
    /** \brief Deserialize into a plugin instance (dispatches on the plugin name) */
    static ProtoFunction* deserialize(DeserializingStream& s);

    /// Documentation string
    static std::string meta_doc;

    /// Plugin factory
    static Function create(const std::string& solver,
                           const std::string& name,
                           const GraphBuilderInternal* gb,
                           const std::vector<std::string>& inputs,
                           const std::vector<std::string>& outputs,
                           const Dict& opts);

    // Infix used to form plugin registration symbols (casadi_register_onnx_<name>)
    static const std::string infix_;

    /// Plugin registry
    static std::map<std::string, Plugin> solvers_;

    /// Mutex protecting the plugin registry
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

   protected:
    /// Reconstruct the frozen metadata from a stream (concrete plugins chain to this)
    explicit OnnxFunction(DeserializingStream& s);

    /// Map an ONNX N-D shape to a 2-D CasADi sparsity (rank<=2 direct, rank>2 flattened column)
    static Sparsity tensor_sparsity(const std::vector<casadi_int>& shape);

    /// Compute the per-model-input feed map: in_src_ (arg index / -2 baked / -1 default) + in_val_
    void build_io_map();

    /// True if input/output index is differentiable (is_diff_in/out, default true)
    bool diff_in(casadi_int i) const { return is_diff_in_.empty() || is_diff_in_.at(i); }
    bool diff_out(casadi_int i) const { return is_diff_out_.empty() || is_diff_out_.at(i); }

    /// Freeze a function from raw model bytes: stage a transient GraphBuilder, then create()
    static Function from_model_data(const std::string& solver, const std::string& name,
                                    const std::vector<uint8_t>& model_data,
                                    const std::vector<std::string>& inputs,
                                    const std::vector<std::string>& outputs,
                                    const Dict& opts);

    /// Build a derivative function with CasADi's full signature by re-create()ing from the model:
    /// present derivative tensors are wired through, absent ones (non-diff) become zero outputs.
    Function wrap_derivative(const std::string& name,
                             const std::vector<std::string>& inames,
                             const std::vector<std::string>& onames,
                             const std::vector<Sparsity>& in_sp,
                             const std::vector<Sparsity>& out_sp,
                             const Dict& dim_bind, const Dict& opts) const;

    /// Serialized ONNX model
    std::vector<uint8_t> model_data_;

    /// Metadata for the exposed inputs/outputs (the selection)
    std::vector<OnnxTensorInfo> in_, out_;

    /// Metadata for every model input (a runtime backend must feed all of them)
    std::vector<OnnxTensorInfo> all_in_;

    /// Per all_in_ entry: exposed-arg index (>=0), -2 baked value, or -1 unwired (default)
    std::vector<casadi_int> in_src_;

    /// Baked input values, flat over all_in_ (numel each; placeholder block when not baked)
    std::vector<double> in_val_;

    /// Names of every input/output tensor in the model (for derivative detection)
    std::set<std::string> model_inputs_, model_outputs_;

    /// Symbolic dimensions naming the forward/adjoint seed counts (bound in get_forward/reverse)
    std::string fwd_dim_ = "nfwd";
    std::string adj_dim_ = "nadj";

    /// Baked input values: input name -> value; such inputs are not exposed as Function inputs
    std::map<std::string, std::vector<double>> input_values_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_ONNX_FUNCTION_IMPL_HPP
