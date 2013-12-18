/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SOLVE_IMPL_HPP
#define SOLVE_IMPL_HPP

#include "solve.hpp"
#include "../fx/linear_solver_internal.hpp"

using namespace std;

namespace CasADi{

  template<bool Tr>
  Solve<Tr>::Solve(const MX& r, const MX& A, const LinearSolver& linear_solver) : linear_solver_(linear_solver){
    casadi_assert_message(r.size2() == A.size1(),"Solve::Solve: dimension mismatch.");
    setDependencies(r,A);
    setSparsity(r.sparsity());
  }

  template<bool Tr>
  void Solve<Tr>::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "(";
    } else if(part==1){
      stream << "/";
    } else {
      if(Tr) stream << "'";
      stream << ")";
    }
  }

  template<bool Tr>
  void Solve<Tr>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    linear_solver_->evaluateDGen(input,output,fwdSeed,fwdSens,adjSeed,adjSens,Tr);
  }

  template<bool Tr>
  void Solve<Tr>::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    linear_solver_->evaluateMXGen(input,output,fwdSeed,fwdSens,adjSeed,adjSens,output_given,Tr);
  }
  
  template<bool Tr>
  void Solve<Tr>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd){
    linear_solver_->propagateSparsityGen(input,output,itmp,rtmp,fwd,Tr);
  }

  template<bool Tr>
  void Solve<Tr>::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    linear_solver_ = deepcopy(linear_solver_, already_copied);
  }

} // namespace CasADi

#endif // SOLVE_IMPL_HPP

