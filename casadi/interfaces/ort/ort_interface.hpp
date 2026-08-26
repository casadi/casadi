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


#ifndef CASADI_ORT_INTERFACE_HPP
#define CASADI_ORT_INTERFACE_HPP

#include "casadi/core/onnx_function_impl.hpp"
#include <casadi/interfaces/ort/casadi_onnx_ort_export.h>

#include <onnxruntime_c_api.h>
#include "ort_runtime.h"

/** \pluginsection{Onnx,ort} */

/// \cond INTERNAL
namespace casadi {

  /** \brief Per-checkout ONNX Runtime state (session + reusable eval scaffolding)

      Holds the mutable, per-evaluation handles so a single OnnxRuntimeInterface can be evaluated
      re-entrantly through independent memory objects, managed by alloc_mem/init_mem/free_mem.

      */
  struct OnnxRuntimeMemory : public OnnxMemory {
    casadi_onnxruntime_data d{};
  };

  /** \brief \pluginbrief{Onnx,ort}

      Black-box evaluation of an ONNX model through Microsoft's ONNX Runtime.

      \author Joris Gillis
      \date 2026

      @copydoc OnnxFunction_doc
      @copydoc plugin_Onnx_ort

      */
  class OnnxRuntimeInterface : public OnnxFunction {
   public:
    OnnxRuntimeInterface(const std::string& name,
                         const GraphBuilderInternal* gb,
                         const std::vector<std::string>& inputs,
                         const std::vector<std::string>& outputs);

    /// Destructor: clear the memory pool HERE (the most-derived class) so the
    /// virtual free_mem dispatches to OnnxRuntimeInterface::free_mem. Relying on
    /// the base ~OnnxFunction::clear_mem() is too late (this subobject is already
    /// destroyed) and would leak every Ort session/env/tensor (idempotent: the
    /// later ~OnnxFunction call then finds an empty pool).
    ~OnnxRuntimeInterface() override;

    /// Plugin factory
    static OnnxFunction* creator(const std::string& name,
                                 const GraphBuilderInternal* gb,
                                 const std::vector<std::string>& inputs,
                                 const std::vector<std::string>& outputs,
                                 const Dict& opts) {
      return new OnnxRuntimeInterface(name, gb, inputs, outputs);
    }

    const char* plugin_name() const override { return "ort"; }
    std::string class_name() const override { return "OnnxRuntimeInterface"; }

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_; }
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new OnnxRuntimeMemory(); }

    /** \brief Initialize memory block: create the ONNX Runtime session */
    int init_mem(void* mem) const override;

    /** \brief Free memory block: tear down the session and reusable scaffolding */
    void free_mem(void* mem) const override;

    /** \brief Evaluate numerically */
    int eval(const double** arg, double** res,
             casadi_int* iw, double* w, void* mem) const override;

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true; }

    /** \brief Generate code for the declarations */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Generate code for the body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into an OnnxRuntimeInterface */
    static ProtoFunction* deserialize(DeserializingStream& s) {
      return new OnnxRuntimeInterface(s);
    }

    /// Documentation
    static const std::string meta_doc;

   private:
    /// Reconstruct from a stream: restore base metadata then rebuild the ORT session
    explicit OnnxRuntimeInterface(DeserializingStream& s);

    /// Abort with the message held by an OrtStatus (and release it)
    void ort_check(OrtStatus* status, const std::string& what) const;

    /// Build the runtime prob_ struct (and its backing arrays) from the frozen base metadata
    void build_prob();

    const OrtApi* ort_api_;
    std::string provider_;

    // Backing storage for prob_ (kept alive for the runtime's lifetime); inputs cover all_in_
    std::vector<const char*> in_names_c_, out_names_c_;
    std::vector<casadi_int> in_elem_, out_elem_, in_ndim_, out_ndim_;
    std::vector<casadi_int> in_dims_, out_dims_, in_numel_, out_numel_;
    casadi_onnxruntime_prob prob_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_ORT_INTERFACE_HPP
