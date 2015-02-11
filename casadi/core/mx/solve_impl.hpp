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
  void Solve<Tr>::evaluateD(const double* const* input, double** output,
                            int* itmp, double* rtmp) {
    // Factorize the matrix
    linear_solver_.setInput(input[1], LINSOL_A);
    linear_solver_.prepare();

    // Solve for nondifferentiated output
    if (input[0]!=input[0]) {
      copy(input[0], input[0]+dep(1).size2(), output[0]);
    }
    linear_solver_.solve(output[0], dep(0).size2(), Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evaluateSX(const SXElement* const* input, SXElement** output,
                             int* itmp, SXElement* rtmp) {
    linear_solver_->evaluateSXGen(input, output, itmp, rtmp, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::eval(const MXPtrV& input, MXPtrV& output) {
    if (input[0]->isZero()) {
      *output[0] = MX(input[0]->shape());
    } else {
      *output[0] = linear_solver_->solve(*input[1], *input[0], Tr);
    }
  }

  template<bool Tr>
  void Solve<Tr>::evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens) {
    linear_solver_->evalFwdLinsol(shared_from_this<MX>(), fwdSeed, fwdSens, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens) {
    linear_solver_->evalAdjLinsol(shared_from_this<MX>(), adjSeed, adjSens, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::propagateSparsity(double** input, double** output,
                                    int* itmp, bvec_t* rtmp, bool fwd) {
    linear_solver_->propagateSparsityGen(input, output, itmp, rtmp, fwd, Tr, dep(0).size2());
  }

  template<bool Tr>
  void Solve<Tr>::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    linear_solver_ = deepcopy(linear_solver_, already_copied);
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SOLVE_IMPL_HPP

