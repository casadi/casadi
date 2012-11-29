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

#include "solve.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include <vector>

using namespace std;

namespace CasADi{

Solve::Solve(const MX& A, const MX& b){
  casadi_warning("MX Solve node is experimental");
  casadi_assert_message(A.size1() == A.size2(),"Linear system must be square");
  casadi_assert_message(A.size1() == b.size1(),"Dimension mismatch");
  setDependencies(A,b);

  // Create the sparsity pattern for the matrix-matrix product
  CRSSparsity spres = sp_dense(A.size1(),b.size2());

  // Save sparsity
  setSparsity(spres);
}

Solve* Solve::clone() const{
  return new Solve(*this);
}

void Solve::printPart(std::ostream &stream, int part) const{
  if(part==0){
    stream << "(";
  } else if(part==1){
    stream << "\\";
  } else {
    stream << ")";
  }
}

void Solve::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  casadi_error("Not implemented");
}

void Solve::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  casadi_error("Not implemented");
}

void Solve::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  casadi_error("Not implemented");
}

void Solve::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
  casadi_error("Not implemented");
}

} // namespace CasADi

