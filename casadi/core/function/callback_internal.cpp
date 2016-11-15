/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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

using namespace std;

namespace casadi {

  #define TRY_CALL(FCN, OBJ, ...) \
  try { \
    casadi_assert_message((OBJ)!=0, "Callback object has been deleted"); \
    return (OBJ)->FCN(__VA_ARGS__);\
  } catch (std::exception& ex) { \
    casadi_error("Error calling \"" CASADI_ASSERT_STR(FCN) "\" for object " \
                 + name_ + ":\n" + std::string(ex.what())); \
  }

  CallbackInternal::
  CallbackInternal(const std::string& name, Callback *self)
    : FunctionInternal(name), self_(self), own_(false) {
  }

  CallbackInternal::~CallbackInternal() {
    if (own_) {
      // Make self a null-pointer and then delete
      // No ownership since reference counter was decreased in Callback::transfer_ownership
      self_->assignNodeNoCount(0);
      delete self_;
    }
  }

  size_t CallbackInternal::get_n_in() {
    TRY_CALL(get_n_in, self_);
  }

  size_t CallbackInternal::get_n_out() {
    TRY_CALL(get_n_out, self_);
  }

  Sparsity CallbackInternal::get_sparsity_in(int i) {
    TRY_CALL(get_sparsity_in, self_, i);
  }

  Sparsity CallbackInternal::get_sparsity_out(int i) {
    TRY_CALL(get_sparsity_out, self_, i);
  }

  std::string CallbackInternal::get_name_in(int i) {
    TRY_CALL(get_name_in, self_, i);
  }

  std::string CallbackInternal::get_name_out(int i) {
    TRY_CALL(get_name_out, self_, i);
  }

  void CallbackInternal::init(const Dict& opts) {
    // Initialize the base classes
    FunctionInternal::init(opts);

    // Initialize this
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    self_->init();
  }

  void CallbackInternal::finalize(const Dict& opts) {
    // Finalize this
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    self_->finalize();

    // Finalize the base classes
    FunctionInternal::finalize(opts);
  }

  void CallbackInternal::
  eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    TRY_CALL(eval, self_, arg, res, iw, w, 0);
  }

  bool CallbackInternal::hasFullJacobian() const {
    TRY_CALL(has_jacobian, self_);
  }

  Function CallbackInternal::
  getFullJacobian(const std::string& name,
                  const std::vector<std::string>& i_names,
                  const std::vector<std::string>& o_names,
                  const Dict& opts) {
    TRY_CALL(get_jacobian, self_, name, opts);
  }

  Function CallbackInternal::
  get_forward(const std::string& name, int nfwd,
              const std::vector<std::string>& i_names,
              const std::vector<std::string>& o_names, const Dict& opts) {
    TRY_CALL(get_forward_new, self_, name, nfwd, i_names, o_names, opts);
  }

  int CallbackInternal::get_n_forward() const {
    TRY_CALL(get_n_forward, self_);
  }

  Function CallbackInternal::
  get_reverse(const std::string& name, int nadj,
              const std::vector<std::string>& i_names,
              const std::vector<std::string>& o_names, const Dict& opts) {
    TRY_CALL(get_reverse_new, self_, name, nadj, i_names, o_names, opts);
  }

  int CallbackInternal::get_n_reverse() const {
    TRY_CALL(get_n_reverse, self_);
  }
} // namespace casadi
