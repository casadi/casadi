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

#include "matrix_matrix_op.hpp"
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

MatrixMatrixOp::MatrixMatrixOp(OPERATION op_, const MX& x, const MX& y) : op(op_){
  setDependencies(x,y);
  casadi_assert_message(x.size1() == y.size1() && x.size2() == y.size2(), "MatrixMatrixOp: dimension mismatch");
  same_sparsity_ = x->sparsity() == y->sparsity();
  casadi_assert_message(same_sparsity_, "different sparsity not implemented");
  setSparsity(x->sparsity());
}

MatrixMatrixOp* MatrixMatrixOp::clone() const{
  return new MatrixMatrixOp(*this);
}

void MatrixMatrixOp::print(std::ostream &stream, const std::vector<std::string>& args) const{
  print_c[op](stream,args.at(0),args.at(1));
}

void MatrixMatrixOp::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  if(same_sparsity_){
    if(nfwd==0 && nadj==0){
      // No sensitivities
      for(int i=0; i<size(); ++i)
        nfun0[op](input[0][i],input[1][i],&output[i]);
      
    } else {
      // Sensitivities
      double tmp[3];  // temporary variable to hold value and partial derivatives of the function
      for(int i=0; i<size(); ++i){
        // Evaluate and get partial derivatives
        nfun1[op](input[0][i],input[1][i],tmp);
        output[i] = tmp[0];
        
        // Propagate forward seeds
        for(int d=0; d<nfwd; ++d){
          fwdSens[d][i] = tmp[1]*fwdSeed[0][d][i] + tmp[2]*fwdSeed[1][d][i];
        }

        // Propagate adjoint seeds
        for(int d=0; d<nadj; ++d){
          adjSens[0][d][i] += adjSeed[d][i]*tmp[1];
          adjSens[1][d][i] += adjSeed[d][i]*tmp[2];
        }
      }
    }
  } else {
    casadi_assert_message(0, "different sparsity not implemented");
  }
}

} // namespace CasADi

