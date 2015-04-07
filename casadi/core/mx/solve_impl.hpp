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
  void Solve<Tr>::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "(";
    } else if (part==1) {
      stream << "'/";
    } else {
      if (!Tr) stream << "'";
      stream << ")'";
    }
  }

  template<bool Tr>
  void Solve<Tr>::print(std::ostream &stream, long& remaining_calls) const {
    if (remaining_calls>0) {
      remaining_calls--;
      stream << "(";
      dep(1)->print(stream, remaining_calls);
      if (Tr) stream << "'";
      stream << "\\";
      dep(0)->print(stream, remaining_calls);
      stream << ")";
    } else {
      stream << "...";
    }
  }

  template<bool Tr>
  void Solve<Tr>::evalD(cp_double* arg, p_double* res,
                        int* itmp, double* rtmp) {
    if (arg[0]!=res[0]) copy(arg[0], arg[0]+dep(0).nnz(), res[0]);
    linear_solver_.setInput(arg[1], LINSOL_A);
    linear_solver_.prepare();
    linear_solver_.solve(res[0], dep(0).size2(), Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evalSX(cp_SXElement* arg, p_SXElement* res, int* itmp, SXElement* rtmp) {
    linear_solver_->evalSXLinsol(arg, res, itmp, rtmp, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    if (arg[0].isZero()) {
      res[0] = MX(arg[0].shape());
    } else {
      res[0] = linear_solver_->solve(arg[1], arg[0], Tr);
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
    linear_solver_->callForwardLinsol(arg, res, fseed, fsens, Tr);
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
    linear_solver_->callReverseLinsol(arg, res, aseed, asens, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::spFwd(cp_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    linear_solver_->spFwdLinsol(arg, res, itmp, rtmp, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::spAdj(p_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    linear_solver_->spAdjLinsol(arg, res, itmp, rtmp, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    linear_solver_ = deepcopy(linear_solver_, already_copied);
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SOLVE_IMPL_HPP

