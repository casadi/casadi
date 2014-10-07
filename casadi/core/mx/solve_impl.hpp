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
  void Solve<Tr>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                            std::vector<double>& rtmp) {
    linear_solver_->evaluateDGen(input, output, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                             std::vector<SXElement>& rtmp) {
    linear_solver_->evaluateSXGen(input, output, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                             MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                             bool output_given) {
    linear_solver_->evaluateMXGen(input, output, fwdSeed, fwdSens, adjSeed, adjSens,
                                  output_given, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                                    std::vector<double>& rtmp, bool fwd) {
    linear_solver_->propagateSparsityGen(input, output, itmp, rtmp, fwd, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    linear_solver_ = deepcopy(linear_solver_, already_copied);
  }

} // namespace casadi

#endif // CASADI_SOLVE_IMPL_HPP

