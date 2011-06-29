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

#include "matrix_scalar_op.hpp"
#include "mx_tools.hpp"
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

MatrixScalarOp::MatrixScalarOp(Operation op_, MX x, const MX& y) : op(op_){
  // Put densifying node in between if necessary
  if(!SX::f0x_is_zero_[op]){
    makeDense(x);
  }
  
  if (y.size()==0) {
    setDependencies(x,0);
  } else {
    setDependencies(x,y);
  }
  setSparsity(x.sparsity());
}

MatrixScalarOp* MatrixScalarOp::clone() const{
  return new MatrixScalarOp(*this);
}

void MatrixScalarOp::print(std::ostream &stream, const std::vector<std::string>& args) const{
  casadi_math<double>::print[op](stream,args.at(0),args.at(1));
}

void MatrixScalarOp::evaluate(const VDptr& input, DMatrix& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  vector<double>& outputd = output.data();

  if(nfwd==0 && nadj==0){
    // No sensitivities
    for(int i=0; i<size(); ++i)
      casadi_math<double>::fun[op](input[0][i],input[1][0],outputd[i]);
    
  } else {
    // Sensitivities
    double tmp[2];  // temporary variable to hold value and partial derivatives of the function
    for(int i=0; i<size(); ++i){
      // Evaluate and get partial derivatives
      casadi_math<double>::fun[op](input[0][i],input[1][0],outputd[i]);
      casadi_math<double>::der[op](input[0][i],input[1][0],outputd[i],tmp);

      // Propagate forward seeds
      for(int d=0; d<nfwd; ++d){
        fwdSens[d][i] = tmp[0]*fwdSeed[0][d][i] + tmp[1]*fwdSeed[1][d][0];
      }

      // Propagate adjoint seeds
      for(int d=0; d<nadj; ++d){
        adjSens[0][d][i] += adjSeed[d][i]*tmp[0];
        adjSens[1][d][0] += adjSeed[d][i]*tmp[1];
      }
    }
  }
}

MX MatrixScalarOp::adFwd(const std::vector<MX>& jx){
  casadi_assert_message(op==SUB || op==ADD || op==MUL, "only addition, subtraction and multiplication implemented (quick hack)");

  if(op==SUB)
    return jx[0]-jx[1];
  else if(op==ADD)
    return jx[0]+jx[1];
  else if(op==MUL)
    return dep(0)*jx[1] + jx[0]*dep(1);
        
  return MX();
}

void MatrixScalarOp::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
  SXMatrix r;
  r.binary(casadi_math<SX>::funE[op],*input[0],*input[1]);
  casadi_assert(output.sparsity()==r.sparsity());
  output.set(r);
}


} // namespace CasADi

