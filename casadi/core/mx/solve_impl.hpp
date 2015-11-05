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
#include "../function/linear_solver_internal.hpp"
/// \cond INTERNAL

using namespace std;

namespace casadi {

  template<bool Tr>
  Solve<Tr>::Solve(const MX& r, const MX& A, const LinearSolver& linear_solver) :
      linear_solver_(linear_solver) {
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
  void Solve<Tr>::evalD(void* mem, const double** arg, double** res,
                        int* iw, double* w) {
    if (arg[0]!=res[0]) copy(arg[0], arg[0]+dep(0).nnz(), res[0]);
    linear_solver_.setInputNZ(arg[1], LINSOL_A);
    linear_solver_.prepare();
    linear_solver_.solve(res[0], dep(0).size2(), Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evalSX(void* mem, const SXElem** arg, SXElem** res, int* iw, SXElem* w) {
    linear_solver_->linsol_evalSX(mem, arg, res, iw, w, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    if (arg[0].is_zero()) {
      res[0] = MX(arg[0].size());
    } else {
      res[0] = linear_solver_->linsol_solve(arg[1], arg[0], Tr);
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
    linear_solver_->linsol_forward(arg, res, fseed, fsens, Tr);
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
    linear_solver_->linsol_reverse(arg, res, aseed, asens, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::spFwd(void* mem, const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    linear_solver_->linsol_spFwd(mem, arg, res, iw, w, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::spAdj(void* mem, bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    linear_solver_->linsol_spAdj(mem, arg, res, iw, w, Tr, dep(0).size2());
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SOLVE_IMPL_HPP

