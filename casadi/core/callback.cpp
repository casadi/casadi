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
    own(new CallbackInternal(name, this));
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

  std::vector<DM> Callback::eval(const std::vector<DM>& arg) const {
    return (*this)->FunctionInternal::eval_dm(arg);
  }

  int Callback::eval(const double** arg, double** res, int* iw, double* w, void* mem) const {
    return (*this)->FunctionInternal::eval(arg, res, iw, w, mem);
  }

  int Callback::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem) const {
    return (*this)->FunctionInternal::eval_sx(arg, res, iw, w, mem);
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

  bool Callback::uses_output() const {
    return (*this)->FunctionInternal::uses_output();
  }

  bool Callback::has_jacobian() const {
    return (*this)->FunctionInternal::has_jacobian();
  }

  Function Callback::
  get_jacobian(const std::string& name,
               const std::vector<std::string>& inames,
               const std::vector<std::string>& onames,
               const Dict& opts) const {
    return (*this)->FunctionInternal::get_jacobian(name, inames, onames, opts);
  }

  Function Callback::
  get_forward(int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    return (*this)->FunctionInternal::get_forward(nfwd, name, inames, onames, opts);
  }

  bool Callback::has_forward(int nfwd) const {
    return (*this)->FunctionInternal::has_forward(nfwd);
  }

  Function Callback::
  get_reverse(int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    return (*this)->FunctionInternal::get_reverse(nadj, name, inames, onames, opts);
  }

  bool Callback::has_reverse(int nadj) const {
    return (*this)->FunctionInternal::has_reverse(nadj);
  }

  void Callback::alloc_w(size_t sz_w, bool persist) {
    return (*this)->alloc_w(sz_w, persist);
  }

  void Callback::alloc_iw(size_t sz_iw, bool persist) {
    return (*this)->alloc_iw(sz_iw, persist);
  }

  void Callback::alloc_arg(size_t sz_arg, bool persist) {
    return (*this)->alloc_arg(sz_arg, persist);
  }

  void Callback::alloc_res(size_t sz_res, bool persist) {
    return (*this)->alloc_res(sz_res, persist);
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
