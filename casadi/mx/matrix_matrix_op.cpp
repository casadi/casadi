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

MatrixMatrixOp::MatrixMatrixOp(Operation op_, const MX& x, const MX& y) : op(op_){
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

void MatrixMatrixOp::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  // Number of non-zeros
  int n = size();
  int nx = dep(0).size();
  int ny = dep(1).size();
      
  if(same_sparsity_){
    if(nfwd==0 && nadj==0){
      // No sensitivities
      for(int i=0; i<n; ++i)
        casadi_math<double>::fun[op](input[0][i],input[1][i],output[i]);
      
    } else {
      // Sensitivities
      double tmp[2];  // temporary variable to hold value and partial derivatives of the function
      for(int i=0; i<n; ++i){
        // Evaluate and get partial derivatives
        casadi_math<double>::fun[op](input[0][i],input[1][i],output[i]);
        casadi_math<double>::der[op](input[0][i],input[1][i],output[i],tmp);
        
        // Propagate forward seeds
        for(int d=0; d<nfwd; ++d){
          fwdSens[d][i] = tmp[0]*fwdSeed[0][d][i] + tmp[1]*fwdSeed[1][d][i];
        }

        // Propagate adjoint seeds
        for(int d=0; d<nadj; ++d){
          adjSens[0][d][i] += adjSeed[d][i]*tmp[0];
          adjSens[1][d][i] += adjSeed[d][i]*tmp[1];
        }
      }
    }
  } else {
    // Non-zero indices of the arguments
    int ix=0, iy=0;
    
    if(nfwd==0 && nadj==0){
      // No sensitivities
      for(int i=0; i<n; ++i){
        double x = mapping_[i]<=0 ? input[0][ix++] : 0;
        double y = mapping_[i]>=0 ? input[1][iy++] : 0;
        casadi_math<double>::fun[op](x,y,output[i]);
      }

    } else {
      // Sensitivities
      double tmp[2];  // temporary variable to hold value and partial derivatives of the function
      for(int i=0; i<n; ++i){
        bool isx = mapping_[i]<=0;
        bool isy = mapping_[i]>=0;

        double x = isx * input[0][ix];
        double y = isy * input[1][iy];

        // Evaluate and get partial derivatives
        casadi_math<double>::fun[op](x,y,output[i]);
        casadi_math<double>::der[op](x,y,output[i],tmp);
        
        // Propagate forward seeds
        for(int d=0; d<nfwd; ++d){
          fwdSens[d][i] = isx*tmp[0]*fwdSeed[0][d][ix] + isy*tmp[1]*fwdSeed[1][d][iy];
        }

        // Propagate adjoint seeds
        for(int d=0; d<nadj; ++d){
          adjSens[0][d][ix] += isx*adjSeed[d][i]*tmp[0];
          adjSens[1][d][iy] += isy*adjSeed[d][i]*tmp[1];
        }
        
        // Increase argument indices
        if(isx && ix+1<nx) ix++; 
        if(isy && iy+1<ny) iy++;
      }

    }
  }
}

MX MatrixMatrixOp::adFwd(const std::vector<MX>& jx){
  casadi_assert_message(op==SUB, "only subtraction implemented (quick hack)");
  return jx[0]-jx[1];
}

} // namespace CasADi

