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

#include "callback_internal.hpp"

namespace casadi {

  #define TRY_CALL(FCN, OBJ, ...) \
  try { \
    casadi_assert((OBJ)!=0, "Callback object has been deleted"); \
    return (OBJ)->FCN(__VA_ARGS__);\
  } catch (std::exception& ex) { \
    casadi_error("Error calling \"" CASADI_STR(FCN) "\" for object " \
                 + name_ + ":\n" + std::string(ex.what())); \
  }

  CallbackInternal::
  CallbackInternal(const std::string& name, Callback *self)
    : FunctionInternal(name), self_(self), has_eval_buffer_(false) {
  }

  CallbackInternal::~CallbackInternal() {
    clear_mem();
  }

  size_t CallbackInternal::get_n_in() {
    TRY_CALL(get_n_in, self_);
  }

  size_t CallbackInternal::get_n_out() {
    TRY_CALL(get_n_out, self_);
  }

  Sparsity CallbackInternal::get_sparsity_in(casadi_int i) {
    TRY_CALL(get_sparsity_in, self_, i);
  }

  Sparsity CallbackInternal::get_sparsity_out(casadi_int i) {
    TRY_CALL(get_sparsity_out, self_, i);
  }

  std::string CallbackInternal::get_name_in(casadi_int i) {
    TRY_CALL(get_name_in, self_, i);
  }

  std::string CallbackInternal::get_name_out(casadi_int i) {
    TRY_CALL(get_name_out, self_, i);
  }

  bool CallbackInternal::has_eval_buffer() const {
    TRY_CALL(has_eval_buffer, self_);
  }

  bool CallbackInternal::has_jac_sparsity(casadi_int oind, casadi_int iind) const {
    TRY_CALL(has_jac_sparsity, self_, oind, iind);
  }

  Sparsity CallbackInternal::get_jac_sparsity(casadi_int oind, casadi_int iind,
      bool symmetric) const {
    TRY_CALL(get_jac_sparsity, self_, oind, iind, symmetric);
  }

  void CallbackInternal::init(const Dict& opts) {
    // Initialize the base classes
    FunctionInternal::init(opts);

    // Initialize this
    casadi_assert(self_!=nullptr, "Callback object has been deleted");
    self_->init();
  }

  void CallbackInternal::finalize() {
    // Finalize this
    casadi_assert(self_!=nullptr, "Callback object has been deleted");
    self_->finalize();

    // Finalize the base classes
    FunctionInternal::finalize();

    has_eval_buffer_ = has_eval_buffer();

    if (has_eval_buffer_) {
      sizes_arg_.resize(n_in_);
      for (casadi_int i=0;i<n_in_;++i) {
        sizes_arg_[i] = nnz_in(i);
      }
      sizes_res_.resize(n_out_);
      for (casadi_int i=0;i<n_out_;++i) {
        sizes_res_[i] = nnz_out(i);
      }
    }
  }

  std::vector<DM> CallbackInternal::eval_dm(const std::vector<DM>& arg) const {
    TRY_CALL(eval, self_, arg);
  }

  /** \brief  Evaluate numerically */
  int CallbackInternal::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    if (has_eval_dm()) {
      return FunctionInternal::eval(arg, res, iw, w, mem);
    } else {
      TRY_CALL(eval_buffer, self_, arg, sizes_arg_, res, sizes_res_);
    }
  }

  bool CallbackInternal::uses_output() const {
    TRY_CALL(uses_output, self_);
  }

  bool CallbackInternal::has_jacobian() const {
    TRY_CALL(has_jacobian, self_);
  }

  Function CallbackInternal::
  get_jacobian(const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
    TRY_CALL(get_jacobian, self_, name, inames, onames, opts);
  }

  Function CallbackInternal::
  get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames, const Dict& opts) const {
    TRY_CALL(get_forward, self_, nfwd, name, inames, onames, opts);
  }

  bool CallbackInternal::has_forward(casadi_int nfwd) const {
    TRY_CALL(has_forward, self_, nfwd);
  }

  Function CallbackInternal::
  get_reverse(casadi_int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames, const Dict& opts) const {
    TRY_CALL(get_reverse, self_, nadj, name, inames, onames, opts);
  }

  bool CallbackInternal::has_reverse(casadi_int nadj) const {
    TRY_CALL(has_reverse, self_, nadj);
  }
} // namespace casadi
