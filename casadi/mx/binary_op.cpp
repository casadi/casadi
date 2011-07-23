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

#include "binary_op.hpp"
#include "mx_tools.hpp"
#include <vector>
#include <sstream>
#include <cassert>
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

BinaryOp::BinaryOp(Operation op, MX x, MX y) : op_(op){
  // If scalar-matrix operation
  if(x.numel()==1 && y.numel()>1){
    x = MX(y.size1(),y.size2(),x);
  }
  
  // If matrix-scalar operation
  if(y.numel()==1 && x.numel()>1){
    y = MX(x.size1(),x.size2(),y);
  }
  
  setDependencies(x,y);
  
  // Get the sparsity pattern
  bool f00_is_zero = casadi_math<double>::f00_is_zero[op_];
  bool f0x_is_zero = casadi_math<double>::f0x_is_zero[op_];
  bool fx0_is_zero = casadi_math<double>::fx0_is_zero[op_];
  CRSSparsity sp = x->sparsity().patternUnion(y->sparsity(),mapping_,f00_is_zero,f0x_is_zero,fx0_is_zero);
  setSparsity(sp);
}

void BinaryOp::print(std::ostream &stream, const std::vector<std::string>& args) const{
  casadi_math<double>::print[op_](stream,args.at(0),args.at(1));
}

void BinaryOp::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  if(nfwd==0 && nadj==0){
    DMatrix::binary_no_alloc(casadi_math<double>::funE[op_],*input[0],*input[1],*output[0],mapping_);
  } else {
    vector<double>& output0 = output[0]->data();
    const vector<double> &input0 = input[0]->data();
    const vector<double> &input1 = input[1]->data();

    // Argument values
    double a[2][2] = {{0,0},{0,0}};

    // Nonzero counters
    int el0=0, el1=0, el=0;

    // Partial derivatives
    double pd[2][2] = {{0,0},{0,0}};
    
    // With sensitivities
    for(int i=0; i<mapping_.size(); ++i){
      // Check which elements are nonzero
      unsigned char m = mapping_[i];
      bool nz0 = m & 1;
      bool nz1 = m & 2;
      bool skip_nz = m & 4;
      
      if(!skip_nz){
      
        // Read the next nonzero 
        if(nz0) a[0][1] = input0[el0];
        if(nz1) a[1][1] = input1[el1];
        
        // Evaluate and get partial derivatives
        casadi_math<double>::fun[op_](a[0][nz0], a[1][nz1],output0[el]);
        casadi_math<double>::der[op_](a[0][nz0], a[1][nz1],output0[el],pd[1]);
        
        // Propagate forward seeds
        for(int d=0; d<nfwd; ++d){
          double s = 0;
          if(nz0) s += pd[nz0][0]*fwdSeed[d][0]->data()[el0];
          if(nz1) s += pd[nz1][1]*fwdSeed[d][1]->data()[el1];
          fwdSens[d][0]->data()[el] = s;
        }
        
        // Propagate adjoint seeds
        for(int d=0; d<nadj; ++d){
          double s = adjSeed[d][0]->data()[el];
          if(nz0) adjSens[d][0]->data()[el0] += s*pd[nz0][0];
          if(nz1) adjSens[d][1]->data()[el1] += s*pd[nz1][1];
        }
        
        // Next nonzero for the output
        el++;
      }
      
      // Go to next nonzero for the arguments
      el0 += nz0;
      el1 += nz1;
    }
  }
}

void BinaryOp::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  Matrix<SX>::binary_no_alloc(casadi_math<SX>::funE[op_],*input[0],*input[1],*output[0],mapping_);
}

void BinaryOp::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens){
  // Evaluate function
  casadi_math<MX>::fun[op_](*input[0],*input[1],*output[0]);

  // Number of forward directions
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  if(nfwd>0 || nadj>0){
    // Get partial derivatives
    MX pd[2];
    casadi_math<MX>::der[op_](*input[0],*input[1],*output[0],pd);
    
    // Chain rule
    for(int d=0; d<nfwd; ++d){
      *fwdSens[d][0] = pd[0]*(*fwdSeed[d][0]) + pd[1]*(*fwdSeed[d][1]);
    }
  }
}

BinaryOp* BinaryOp::clone() const{
  return new BinaryOp(*this);
}

} // namespace CasADi

