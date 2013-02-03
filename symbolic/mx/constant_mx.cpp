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

#include "constant_mx.hpp"
#include <cassert>
#include <vector>
#include <algorithm>
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

ConstantMX::ConstantMX(const Matrix<double> &x) : x_(x){
  setSparsity(x.sparsity());
}

ConstantMX* ConstantMX::clone() const{
  return new ConstantMX(*this);
}

void ConstantMX::printPart(std::ostream &stream, int part) const{
  x_.print(stream);
}

void ConstantMX::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  output[0]->set(x_);
  for(int d=0; d<nfwd; ++d){
    fwdSens[d][0]->setZero();
  }
}

bool ConstantMX::__nonzero__() const {
  if (numel()!=1) casadi_error("Can only determine truth value of scalar MX.");
  if (size()!=1) casadi_error("Can only determine truth value of dense scalar MX.");
  return x_.at(0)!=0;
}

void ConstantMX::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  SXMatrix r(x_);
  casadi_assert(output[0]->sparsity()==r.sparsity());
  output[0]->set(r);
}

void ConstantMX::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){

  // Evaluate nondifferentiated
  if(!output_given){
    if(output[0]){
      *output[0] = x_;
    }
  }
  
  // Number of derivative directions
  int nfwd = fwdSens.size();
  if(nfwd==0) return; // Quick return
  
  // Derivatives
  MX zero_sens = MX::sparse(size1(),size2());
  for(int d=0; d<nfwd; ++d){
    if(fwdSens[d][0]){
      *fwdSens[d][0] = zero_sens;
    }
  }
}

void ConstantMX::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
  if(fwd){
    bvec_t *outputd = get_bvec_t(output[0]->data());
    fill_n(outputd,output[0]->size(),0);
  }
}

void ConstantMX::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
  // Print all nonzeros row by row
  for(int k=0; k<x_.size(); ++k){
    
    // Print left hand side of assignment
    stream << "  " << res.front() << "[" << k << "]=";
        
    // Print right hand side of assignment
    stream << x_.at(k) << ";" << endl;
  }
}


} // namespace CasADi

