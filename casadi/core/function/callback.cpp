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

  Callback::Callback() {
  }

  Callback::Callback(const Callback& obj) : Function(obj) {
    casadi_error("Callback objects cannot be copied");
  }

  void Callback::construct(const std::string& name, const Dict& opts) {
    assignNode(new CallbackInternal(name, this));
    (*this)->construct(opts);
  }

  Callback::~Callback() {
    // Make sure that this object isn't used after its deletion
    if (!is_null()) {
      (*this)->self_ = 0;
    }
  }

  Function Callback::create(const std::string& name, Callback* n, const Dict& opts) {
    n->construct(name, opts);
    Function ret = *n;
    n->transfer_ownership();
    return ret;
  }

  std::vector<DM> Callback::eval(const std::vector<DM>& arg) {
    casadi_error("Callback::eval has not been implemented");
    return std::vector<DM>();
  }

  void Callback::eval(const double** arg, double** res, int* iw, double* w, int mem) {
    // Allocate input matrices
    int n_in = this->n_in();
    std::vector<DM> argv(n_in);
    for (int i=0; i<n_in; ++i) {
      argv[i] = DM(sparsity_in(i));
      casadi_copy(arg[i], argv[i].nnz(), argv[i].ptr());
    }

    // Evaluate
    std::vector<DM> resv = eval(argv);

    casadi_assert_message(resv.size()==n_out(),
      "Callback::eval: expected " + to_string(n_out()) + " outputs, got "
      + to_string(resv.size()) +".");

    for (int i=0; i<n_out(); ++i) {
      casadi_assert_message(resv[i].sparsity()==sparsity_out(i),
        "Callback::eval: Shape mismatch for output " << i << ": got " + resv[i].dim() +
        ", expected " + sparsity_out(i).dim() + ".");
    }

    // Get the outputs
    int n_out = this->n_out();
    for (int i=0; i<n_out; ++i) {
      casadi_copy(resv[i].ptr(), resv[i].nnz(), res[i]);
    }
  }

  const CallbackInternal* Callback::operator->() const {
    return static_cast<const CallbackInternal*>(Function::operator->());
  }

  CallbackInternal* Callback::operator->() {
    return static_cast<CallbackInternal*>(Function::operator->());
  }

  int Callback::get_n_in() {
    return (*this)->FunctionInternal::get_n_in();
  }

  int Callback::get_n_out() {
    return (*this)->FunctionInternal::get_n_out();
  }

  Sparsity Callback::get_sparsity_in(int i) {
    return (*this)->FunctionInternal::get_sparsity_in(i);
  }

  Sparsity Callback::get_sparsity_out(int i) {
    return (*this)->FunctionInternal::get_sparsity_out(i);
  }

  std::string Callback::get_name_in(int i) {
    return (*this)->FunctionInternal::get_name_in(i);
  }

  std::string Callback::get_name_out(int i) {
    return (*this)->FunctionInternal::get_name_out(i);
  }

  bool Callback::has_jacobian() const {
    return (*this)->FunctionInternal::hasFullJacobian();
  }

  Function Callback::get_jacobian(const std::string& name, const Dict& opts) {
    return (*this)->FunctionInternal::getFullJacobian(name, opts);
  }

  Function Callback::get_forward(const std::string& name, int nfwd, Dict& opts) {
    return (*this)->FunctionInternal::get_forward_old(name, nfwd, opts);
  }

  Function Callback::get_forward_new(const std::string& name, int nfwd, Dict& opts) {
    return (*this)->FunctionInternal::get_forward(name, nfwd, opts);
  }

  int Callback::get_n_forward() const {
    return (*this)->FunctionInternal::get_n_forward();
  }

  Function Callback::get_reverse(const std::string& name, int nadj, Dict& opts) {
    return (*this)->FunctionInternal::get_reverse_old(name, nadj, opts);
  }

  Function Callback::get_reverse_new(const std::string& name, int nadj, Dict& opts) {
    return (*this)->FunctionInternal::get_reverse(name, nadj, opts);
  }

  int Callback::get_n_reverse() const {
    return (*this)->FunctionInternal::get_n_reverse();
  }

  void Callback::transfer_ownership() {
    casadi_assert_message(!is_null(), "Null pointer.");
    casadi_assert_message(!(*this)->own_, "Ownership has already been transferred.");
    casadi_assert_message(getCount()>1, "There are no owning references");
    // Decrease the reference counter to offset the effect of the owning reference
    count_down();
    // Mark internal class as owning
    (*this)->own_ = true;
  }

} // namespace casadi
