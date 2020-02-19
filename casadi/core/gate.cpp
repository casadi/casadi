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


#include "gate.hpp"
using namespace std;
namespace casadi {


  static std::vector<MX> casadi_gates;

  std::vector<MX>& Gate::gates() {
    return casadi_gates;
  }

  void Gate::gates_clean() {
    casadi_gates.clear();
  }

  
  MX Gate::create(const MX& dep) {
    if (dep.is_gate()) {
      return dep;
    } else {
      MX ret = MX::create(new Gate(dep));
      gates().push_back(ret);
//      std::swap(dep.get(), ret.get());
//      std::swap(dep, ret);
      return ret;
    }
  }

  Gate::Gate(const MX& x) {
    set_dep(x);
    set_sparsity(x.sparsity());
    value_.resize(x.nnz(), nan);
  }

  std::string Gate::disp(const std::vector<std::string>& arg) const {
    return "gate()";
  }

  void Gate::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_gate();
  }

  void Gate::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    MX zero_sens(size1(), size2());
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = zero_sens;
    }

    return;
  }

  void Gate::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    return;
  }

  int Gate::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg) {
        casadi_copy(arg[0], nnz(), get_ptr(value_));
    }
    if (res[0]) casadi_copy(get_ptr(value_), nnz(), res[0]);
    return 0;
  }

  int Gate::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return 0;
  }

  int Gate::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  int Gate::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  void Gate::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("Gate::value", value_);
  }

  Gate::Gate(DeserializingStream& s) : MXNode(s) {
    s.unpack("Gate::value", value_);
  }

  void Gate::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const {
    if (arg.empty()) {
      g << g.copy(g.rw_double(this), nnz(), g.work(res[0], nnz())) << '\n';
    } else {
      g << g.copy(g.work(arg[0], nnz()), nnz(), g.rw_double(this)) << '\n';
      g << g.copy(g.rw_double(this), nnz(), g.work(res[0], nnz())) << '\n';      
    }
  }

  void Gate::add_dependency(CodeGenerator& g) const {
    g.define_rw_double(this, nnz());
  }

  const std::vector<double>& Gate::get_double_vec() const {
    return value_;
  }

} // namespace casadi
