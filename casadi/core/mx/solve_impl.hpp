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


#ifndef CASADI_SOLVE_IMPL_HPP
#define CASADI_SOLVE_IMPL_HPP

#include "solve.hpp"
#include "../function/function_internal.hpp"

using namespace std;

namespace casadi {

  template<bool Tr>
  Solve<Tr>::Solve(const MX& r, const MX& A, const Function& linear_solver) :
      linsol_(linear_solver) {
    casadi_assert_message(r.size1() == A.size2(), "Solve::Solve: dimension mismatch.");
    setDependencies(r, A);
    setSparsity(r.sparsity());
  }

  template<bool Tr>
  std::string Solve<Tr>::print(const std::vector<std::string>& arg) const {
    std::stringstream ss;
    ss << "(" << arg.at(1);
    if (Tr) ss << "'";
    ss << "\\" << arg.at(0) << ")";
    return ss.str();
  }

  template<bool Tr>
  void Solve<Tr>::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    if (arg[0]!=res[0]) copy(arg[0], arg[0]+dep(0).nnz(), res[0]);
    linsol_.setup(arg+2, res+1, iw, w, mem);
    linsol_.linsol_factorize(arg[1], mem);
    linsol_.linsol_solve(res[0], dep(0).size2(), Tr, mem);
  }

  template<bool Tr>
  void Solve<Tr>::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    linsol_->linsol_eval_sx(arg, res, iw, w, mem, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    if (arg[0].is_zero()) {
      res[0] = MX(arg[0].size());
    } else {
      res[0] = linsol_->linsol_solve(arg[1], arg[0], Tr);
    }
  }

  template<bool Tr>
  void Solve<Tr>::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Call the cached functions
    linsol_->linsol_forward(arg, res, fseed, fsens, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Call the cached functions
    linsol_->linsol_reverse(arg, res, aseed, asens, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    linsol_->linsol_spFwd(arg, res, iw, w, mem, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    linsol_->linsol_spAdj(arg, res, iw, w, mem, Tr, dep(0).size2());
  }

  template<bool Tr>
  size_t Solve<Tr>::sz_arg() const {
    return ndep() + linsol_.sz_arg();
  }

  template<bool Tr>
  size_t Solve<Tr>::sz_res() const {
    return nout() + linsol_.sz_res();
  }

  template<bool Tr>
  size_t Solve<Tr>::sz_iw() const {
    return linsol_.sz_iw();
  }

  template<bool Tr>
  size_t Solve<Tr>::sz_w() const {
    return linsol_.sz_w() + sparsity().size1();
  }

} // namespace casadi

#endif // CASADI_SOLVE_IMPL_HPP
