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

#ifndef CASADI_UNO_INTERFACE_HPP
#define CASADI_UNO_INTERFACE_HPP

#include <casadi/interfaces/uno/casadi_nlpsol_uno_export.h>
#include "casadi/core/nlpsol_impl.hpp"
#include <Uno_C_API.h>

namespace casadi {
  #include "uno_runtime.hpp"
}

/** \defgroup plugin_Nlpsol_uno Title
    \par
    David Kiessling
    Uno interface

    \identifier{2j2} */

/** \pluginsection{Nlpsol,uno} */

/// \cond INTERNAL
namespace casadi {
  class UnoInterface;

  struct CASADI_NLPSOL_UNO_EXPORT UnoMemory : public NlpsolMemory {
    const UnoInterface& self;
    casadi_uno_data<double> d_uno;
    const char* return_status;
    UnoMemory(const UnoInterface& uno_interface);
    ~UnoMemory();
  };

  /** \brief \pluginbrief{Nlpsol,uno}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_uno
  */
  class CASADI_NLPSOL_UNO_EXPORT UnoInterface : public Nlpsol {
  public:
    explicit UnoInterface(const std::string& name, const Function& nlp);
    ~UnoInterface() override;

    const char* plugin_name() const override { return "uno"; }
    std::string class_name() const override { return "UnoInterface"; }

    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new UnoInterface(name, nlp);
    }

    static const Options options_;
    const Options& get_options() const override { return options_; }

    void init(const Dict& opts) override;
    void* alloc_mem() const override { return new UnoMemory(*this); }
    int  init_mem(void* mem) const override;
    void free_mem(void* mem) const override;
    void set_work(void* mem, const double**& arg, double**& res,
                  casadi_int*& iw, double*& w) const override;
    int  solve(void* mem) const override;
    Dict get_stats(void* mem) const override;

    void serialize_body(SerializingStream& s) const override;
    static ProtoFunction* deserialize(DeserializingStream& s) {
      return new UnoInterface(s);
    }

    static const std::string meta_doc;

    // Codegen
    void codegen_body(CodeGenerator& g) const override;
    void codegen_declarations(CodeGenerator& g) const override;
    void codegen_init_mem(CodeGenerator& g) const override;
    void codegen_free_mem(CodeGenerator& g) const override;
    std::string codegen_mem_type() const override { return "struct casadi_uno_data"; }
    bool codegen_needs_mem() const override { return true; }

  protected:
    explicit UnoInterface(DeserializingStream& s);

  private:
    // Build the prob struct (sparsity ptrs + callback fn ptrs) from members.
    // Used at C++ init() time. Codegen has set_uno_prob(g) instead.
    void set_uno_prob();
    void set_uno_prob(CodeGenerator& g) const;

    // NLP sparsities discovered in init().
    Sparsity jacg_sp_;
    Sparsity hesslag_sp_;
    // Sparsity index arrays in Uno's int width, filled once in init().
    std::vector<uno_int> jacobian_row_indices_;
    std::vector<uno_int> jacobian_column_indices_;
    std::vector<uno_int> hessian_row_indices_;
    std::vector<uno_int> hessian_column_indices_;
    // -inf / +inf placeholder constraint bounds passed to uno_set_constraints at
    // init_mem time -- kept as members so the pointer stays valid across the call
    // (real bounds are applied per-solve). Variable bounds aren't part of the
    // (unconstrained) model, so no x placeholders are needed.
    std::vector<double> placeholder_lb_g_, placeholder_ub_g_;
    // Solver-specific options forwarded to uno (the {"uno": {...}} dict).
    Dict opts_;
    // Filled in by set_uno_prob(). One member per UnoInterface; the prob
    // pointers (sparsities, fn ptrs) outlive every solve.
    casadi_uno_prob<double> p_uno_;
  };

}  // namespace casadi
/// \endcond
#endif  // CASADI_UNO_INTERFACE_HPP
