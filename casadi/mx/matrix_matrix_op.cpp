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
#include "mx_tools.hpp"
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

MatrixMatrixOp::MatrixMatrixOp(Operation op_, MX x, MX y) : op(op_){
  // Put densifying nodes in between if necessary
  if(!SX::f00_is_zero_[op] && !(x.dense() || y.dense())){
    if(x.size()>y.size()) // enough to make one of them dense
      makeDense(x);
    else
      makeDense(y);
  }
  
  setDependencies(x,y);
  
  // Check if the sparsity patterns are the same
  same_sparsity_ = x.sparsity() == y.sparsity();
  
  // Get the sparsity pattern
  if(same_sparsity_){
    setSparsity(x->sparsity());
  } else {
    CRSSparsity sp = x->sparsity().patternUnion(y->sparsity(),mapping_);
    setSparsity(sp);
  }
}

MatrixMatrixOp* MatrixMatrixOp::clone() const{
  return new MatrixMatrixOp(*this);
}

void MatrixMatrixOp::print(std::ostream &stream, const std::vector<std::string>& args) const{
  casadi_math<double>::print[op](stream,args.at(0),args.at(1));
}

void MatrixMatrixOp::evaluate(const std::vector<DMatrix*>& input, DMatrix& output, const vvDMatrixP& fwdSeed, std::vector<DMatrix*>& fwdSens, const std::vector<DMatrix*>& adjSeed, vvDMatrixP& adjSens, int nfwd, int nadj){
  vector<double>& outputd = output.data();
  const vector<double> &input0 = input[0]->data();
  const vector<double> &input1 = input[1]->data();

  // Number of non-zeros
  int n = size();
  int nx = dep(0).size();
  int ny = dep(1).size();
      
  if(same_sparsity_){
    if(nfwd==0 && nadj==0){
      // No sensitivities
      for(int i=0; i<n; ++i)
        casadi_math<double>::fun[op](input0[i],input1[i],outputd[i]);
      
    } else {
      // Sensitivities
      double tmp[2];  // temporary variable to hold value and partial derivatives of the function
      for(int i=0; i<n; ++i){
        // Evaluate and get partial derivatives
        casadi_math<double>::fun[op](input0[i],input1[i],outputd[i]);
        casadi_math<double>::der[op](input0[i],input1[i],outputd[i],tmp);
        
        // Propagate forward seeds
        for(int d=0; d<nfwd; ++d){
          fwdSens[d]->data()[i] = tmp[0]*fwdSeed[0][d]->data()[i] + tmp[1]*fwdSeed[1][d]->data()[i];
        }

        // Propagate adjoint seeds
        for(int d=0; d<nadj; ++d){
          adjSens[0][d]->data()[i] += adjSeed[d]->data()[i]*tmp[0];
          adjSens[1][d]->data()[i] += adjSeed[d]->data()[i]*tmp[1];
        }
      }
    }
  } else {
    // Non-zero indices of the arguments
    int ix=0, iy=0;
    
    if(nfwd==0 && nadj==0){
      // No sensitivities
      for(int i=0; i<n; ++i){
        double x = mapping_[i]<=0 ? input0[ix++] : 0;
        double y = mapping_[i]>=0 ? input1[iy++] : 0;
        casadi_math<double>::fun[op](x,y,outputd[i]);
      }

    } else {
      // Sensitivities
      double tmp[2];  // temporary variable to hold value and partial derivatives of the function
      for(int i=0; i<n; ++i){
        bool isx = mapping_[i]<=0;
        bool isy = mapping_[i]>=0;

        double x = isx ? input0[ix] : 0;
        double y = isy ? input1[iy] : 0;

        // Evaluate and get partial derivatives
        casadi_math<double>::fun[op](x,y,outputd[i]);
        casadi_math<double>::der[op](x,y,outputd[i],tmp);
        
        // Propagate forward seeds
        for(int d=0; d<nfwd; ++d){
          fwdSens[d]->data()[i] = (isx ? tmp[0]*fwdSeed[0][d]->data()[ix]: 0) + (isy ? tmp[1]*fwdSeed[1][d]->data()[iy]: 0);
        }

        // Propagate adjoint seeds
        for(int d=0; d<nadj; ++d){
          if (isx) adjSens[0][d]->data()[ix] += adjSeed[d]->data()[i]*tmp[0];
          if (isy) adjSens[1][d]->data()[iy] += adjSeed[d]->data()[i]*tmp[1];
        }
        
        // Increase argument indices
        if(isx && ix+1<nx) ix++; 
        if(isy && iy+1<ny) iy++;
      }

    }
  }
}

MX MatrixMatrixOp::adFwd(const std::vector<MX>& jx){
  casadi_assert_message(op==SUB || op==ADD || op==MUL, "only addition, subtraction and multiplication implemented (quick hack)");

  if(op==SUB)
    return jx[0]-jx[1];
  else if(op==ADD)
    return jx[0]+jx[1];
  else if(op==MUL)
    return dep(0)*jx[1] + jx[0]*dep(1);
        
  return MX();
}

void MatrixMatrixOp::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
  SXMatrix r;
  r.binary_old(casadi_math<SX>::funE[op],*input[0],*input[1]);
  casadi_assert(output.sparsity()==r.sparsity());
  output.set(r);
}


} // namespace CasADi

