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


#ifndef CASADI_EXTERNAL_IMPL_HPP
#define CASADI_EXTERNAL_IMPL_HPP

#include "external.hpp"
#include "function_internal.hpp"

#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#endif
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL

/// \cond INTERNAL

namespace casadi {

class CASADI_EXPORT External : public FunctionInternal {
 protected:
  /** \brief Information about the library

      \identifier{1z2} */
  Importer li_;

  /** \brief Increase/decrease reference counter

      \identifier{1z3} */
  signal_user_data_t incref_, decref_;

  /** \brief Number of inputs and outputs

      \identifier{1z4} */
  getint_user_data_t get_n_in_, get_n_out_;

  /** \brief Names of inputs and outputs

      \identifier{1z5} */
  name_user_data_t get_name_in_, get_name_out_;

  /** \brief Get default inputs

      \identifier{1z6} */
  default_user_data_t get_default_in_;

  /** \brief Work vector sizes

      \identifier{1z7} */
  work_user_data_t work_;

  ///@{
  /** \brief Data vectors

      \identifier{1z8} */
  std::vector<casadi_int> int_data_;
  std::vector<double> real_data_;
  std::string string_data_;
  ///@}
 public:

  /** \brief Constructor

      \identifier{1z9} */
  External(const std::string& name, const Importer& li);

  /** \brief Destructor

      \identifier{1za} */
  ~External() override = 0;

  /** \brief Any symbol found?

      \identifier{1zb} */
  virtual bool any_symbol_found() const;

  // Factory
  Function factory(const std::string& name,
                           const std::vector<std::string>& s_in,
                           const std::vector<std::string>& s_out,
                           const Function::AuxOut& aux,
                           const Dict& opts) const override;

  /** \brief Get type name

      \identifier{1zc} */
  std::string class_name() const override { return "External";}

  /** \brief Initialize members that are unique

      \identifier{1zd} */
  virtual void init_external();

  /// Initialize
  void init(const Dict& opts) override;

  void init_post_options() override;

  /** \brief Generate code for the declarations of the C function

      \identifier{1ze} */
  void codegen_declarations(CodeGenerator& g) const override;

  /** \brief Generate code for the body of the C function

      \identifier{1zf} */
  void codegen_body(CodeGenerator& g) const override;

  /** \brief Thread-local memory object type

      \identifier{1zg} */
  std::string codegen_mem_type() const override;

  /** \brief Codegen incref for dependencies

      \identifier{1zh} */
  void codegen_incref(CodeGenerator& g) const override;

  /** \brief Codegen decref for dependencies

      \identifier{1zi} */
  void codegen_decref(CodeGenerator& g) const override;

  /** \brief Codegen for checkout

      \identifier{1zj} */
  void codegen_checkout(CodeGenerator& g) const override;

  /** \brief Codegen for release

      \identifier{1zk} */
  void codegen_release(CodeGenerator& g) const override;

  /** \brief Codegen decref for alloc_mem

      \identifier{1zl} */
  void codegen_alloc_mem(CodeGenerator& g) const override;

  /** \brief Codegen decref for init_mem

      \identifier{1zm} */
  void codegen_init_mem(CodeGenerator& g) const override;

  /** \brief Codegen for free_mem

      \identifier{1zn} */
  void codegen_free_mem(CodeGenerator& g) const override;

  ///@{
  /** \brief Number of function inputs and outputs

      \identifier{1zo} */
  size_t get_n_in() override;
  size_t get_n_out() override;
  ///@}

  /** \brief Default inputs

      \identifier{1zp} */
  double get_default_in(casadi_int i) const override;

  ///@{
  /** \brief Names of function input and outputs

      \identifier{1zq} */
  std::string get_name_in(casadi_int i) override;
  std::string get_name_out(casadi_int i) override;
  /// @}

  ///@{
  /** \brief Forward mode derivatives

      \identifier{1zr} */
  Function get_forward(casadi_int nfwd, const std::string& name,
                               const std::vector<std::string>& inames,
                               const std::vector<std::string>& onames,
                               const Dict& opts) const override;
  bool has_forward(casadi_int nfwd) const override;
  ///@}

  ///@{
  /** \brief Reverse mode derivatives

      \identifier{1zs} */
  Function get_reverse(casadi_int nadj, const std::string& name,
                               const std::vector<std::string>& inames,
                               const std::vector<std::string>& onames,
                               const Dict& opts) const override;
  bool has_reverse(casadi_int nadj) const override;
  ///@}

  ///@{
  /** \brief Full Jacobian

      \identifier{1zt} */
  bool has_jacobian() const override;
  Function get_jacobian(const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const override;
  ///@}

  /** \brief Serialize an object without type information

      \identifier{1zu} */
  void serialize_body(SerializingStream &s) const override;

  /** \brief Deserialize into MX

      \identifier{1zv} */
  static ProtoFunction* deserialize(DeserializingStream& s);

  /** \brief String used to identify the immediate FunctionInternal subclass

      \identifier{1zw} */
  std::string serialize_base_function() const override { return "External"; }

 protected:
  /** \brief Deserializing constructor

      \identifier{1zx} */
  explicit External(DeserializingStream& s);
};

class CASADI_EXPORT GenericExternal : public External {
  // Sparsities
  sparsity_user_data_t get_sparsity_in_, get_sparsity_out_, get_jac_sparsity_;
  // Differentiability
  diff_user_data_t get_diff_in_, get_diff_out_;

 public:
  /** \brief Constructor

      \identifier{1zy} */
  GenericExternal(const std::string& name, const Importer& li);

  /** \brief  Destructor

      \identifier{1zz} */
  ~GenericExternal() override { this->clear_mem();}

  /** \brief Any symbol found?

      \identifier{200} */
  bool any_symbol_found() const override;

  /** \brief Initialize members that are unique

      \identifier{201} */
  void init_external() override;

  /// Initialize
  void init(const Dict& opts) override;

  /// @{
  /** \brief Retreive sparsities

      \identifier{202} */
  Sparsity get_sparsity_in(casadi_int i) override;
  Sparsity get_sparsity_out(casadi_int i) override;
  /// @}

  ///@{
  /** \brief Return sparsity of Jacobian of an output respect to an input

      \identifier{203} */
  bool has_jac_sparsity(casadi_int oind, casadi_int iind) const override;
  Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override;
  ///@}

  /// @{
  /** \brief Retreive differentiability

      \identifier{204} */
  bool get_diff_in(casadi_int i) override;
  bool get_diff_out(casadi_int i) override;
  /// @}

  /** \brief Serialize type information

      \identifier{205} */
  void serialize_type(SerializingStream &s) const override;

  /** \brief Deserializing constructor

      \identifier{206} */
  explicit GenericExternal(DeserializingStream& s);

};

} // namespace casadi
/// \endcond

#endif // CASADI_EXTERNAL_IMPL_HPP
