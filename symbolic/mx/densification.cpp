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

#include "densification.hpp"
#include "mx_tools.hpp"
#include <vector>
#include <sstream>
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

Densification::Densification(const MX& x){
  setDependencies(x);
  setSparsity(CRSSparsity(x.size1(),x.size2(),true));
}

Densification* Densification::clone() const{
  return new Densification(*this);
}

void Densification::printPart(std::ostream &stream, int part) const{
  if(part==0){
    stream << "dense(";
  } else {
    stream << ")";
  }
}

void Densification::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();

  // Propate values
  input[0]->get(output[0]->data(),DENSE);
  
  // Propagate forward seeds
  for(int d=0; d<nfwd; ++d){
    fwdSeed[d][0]->get(fwdSens[d][0]->data(),DENSE);
  }

  // Propagate adjoint seeds
  for(int d=0; d<nadj; ++d){
    adjSens[d][0]->set(adjSeed[d][0]->data(),DENSE);
    adjSeed[d][0]->setZero();
  }
}

void Densification::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  input[0]->get(output[0]->data(),DENSE);
}

void Densification::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  // Evaluate function
  if(!output_given){
    *output[0] = *input[0];
    makeDense(*output[0]);
  }
  
  // Propagate forward seeds
  int nfwd = fwdSens.size();
  for(int d=0; d<nfwd; ++d){
    *fwdSens[d][0] = *fwdSeed[d][0];
  }

  // Propagate adjoint seeds
  int nadj = adjSeed.size();
  for(int d=0; d<nadj; ++d){
    *adjSens[d][0] += *adjSeed[d][0];
    *adjSeed[d][0] = MX();
  }
}

void Densification::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
  bvec_t *inputd = get_bvec_t(input[0]->data());
  bvec_t *outputd = get_bvec_t(output[0]->data());
  
  int d1 = input[0]->size1();
  int d2 = input[0]->size2();
  const vector<int>& rowind = input[0]->rowind();
  const vector<int>& col = input[0]->col();
  
  int k=0; // index of the result
  for(int i=0; i<d1; ++i){ // loop over rows
    for(int el=rowind[i]; el<rowind[i+1]; ++el){ // loop over the non-zero elements
      int j=col[el];  // column
      for(; k<i*d2+j; ++k){
	outputd[k] = 0; // add zeros before the non-zero element
      }
      
      // add the non-zero element
      if(fwd){
        outputd[k] = inputd[el];
      } else {
        inputd[el] |= outputd[k];
	outputd[k] = 0;
      }
      k++;
    }
  }
  
  // add sparse zeros at the end of the matrix
  for(; k<d1*d2; ++k)
    outputd[k] = 0;
}


} // namespace CasADi

