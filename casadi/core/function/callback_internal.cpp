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

  size_t CallbackInternal::get_n_in() const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_n_in();
  }

  size_t CallbackInternal::get_n_out() const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_n_out();
  }

  Sparsity CallbackInternal::get_sparsity_in(int ind) const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_input_sparsity(ind);
  }

  Sparsity CallbackInternal::get_sparsity_out(int ind) const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_output_sparsity(ind);
  }

  void CallbackInternal::init(const Dict& opts) {
    // Initialize the base classes
    FunctionInternal::init(opts);

    // Initialize this
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    self_->init();
  }

  void CallbackInternal::finalize() {
    // Finalize this
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    self_->finalize();

    // Finalize the base classes
    FunctionInternal::finalize();
  }

  void CallbackInternal::
  eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    self_->eval(arg, res, iw, w, 0);
  }

  bool CallbackInternal::hasFullJacobian() const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->has_jacobian();
  }

  Function CallbackInternal::getFullJacobian(const std::string& name, const Dict& opts) {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_jacobian(name, opts);
  }

  Function CallbackInternal::get_forward(const std::string& name, int nfwd, Dict& opts) {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_forward(name, nfwd, opts);
  }

  int CallbackInternal::get_n_forward() const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_n_forward();
  }

  Function CallbackInternal::get_reverse(const std::string& name, int nadj, Dict& opts) {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_reverse(name, nadj, opts);
  }

  int CallbackInternal::get_n_reverse() const {
    casadi_assert_message(self_!=0, "Callback object has been deleted");
    return self_->get_n_reverse();
  }
} // namespace casadi

