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

#include "multiplication.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include <vector>

using namespace std;

namespace CasADi{

Multiplication::Multiplication(const MX& x, const MX& y_trans){
  casadi_assert_message(x.size2() == y_trans.size2(),"Multiplication::Multiplication: dimension mismatch");
  setDependencies(x,y_trans);

  // Create the sparsity pattern for the matrix-matrix product
  CRSSparsity spres = x->sparsity().patternProduct(y_trans.sparsity());

  // Save sparsity
  setSparsity(spres);
}

Multiplication* Multiplication::clone() const{
  return new Multiplication(*this);
}

void Multiplication::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "prod(" << args.at(0) << "," << args.at(1) << ")";
}

void Multiplication::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();

  fill(output[0]->begin(),output[0]->end(),0);
  DMatrix::prod_no_alloc(*input[0],*input[1],*output[0]);

  // Forward sensitivities: dot(Z) = dot(X)*Y + X*dot(Y)
  for(int d=0; d<nfwd; ++d){
    fill(fwdSens[d][0]->begin(),fwdSens[d][0]->end(),0);
    DMatrix::prod_no_alloc(*fwdSeed[d][0],*input[1],*fwdSens[d][0]);
    DMatrix::prod_no_alloc(*input[0],*fwdSeed[d][1],*fwdSens[d][0]);
  }

  // Adjoint sensitivities
  for(int d=0; d<nadj; ++d){
    DMatrix::prod_no_alloc1(*adjSens[d][0],*input[1],*adjSeed[d][0]);
    DMatrix::prod_no_alloc2(*input[0],*adjSens[d][1],*adjSeed[d][0]);
  }
}

void Multiplication::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  fill(output[0]->begin(),output[0]->end(),0);
  SXMatrix::prod_no_alloc(*input[0],*input[1],*output[0]);
}

void Multiplication::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  if(!output_given)
    *output[0] = prod(*input[0],trans(*input[1]));
  int nfwd = fwdSens.size();
  for(int d=0; d<nfwd; ++d){
    *fwdSens[d][0] = prod(*fwdSeed[d][0],trans(*input[1])) + prod(*input[0],trans(*fwdSeed[d][1]));
  }
}

void Multiplication::propagateSparsity(const DMatrixPtrV& input, DMatrixPtrV& output){
  const bvec_t *x_data = get_bvec_t(input[0]->data());
  const bvec_t *y_trans_data = get_bvec_t(input[1]->data());
  bvec_t *z_data = get_bvec_t(output[0]->data());
  
  // Direct access to the arrays
  const std::vector<int> &z_col = output[0]->col();
  const std::vector<int> &z_rowind = output[0]->rowind();
  const std::vector<int> &x_col = input[0]->col();
  const std::vector<int> &y_row = input[1]->col();
  const std::vector<int> &x_rowind = input[0]->rowind();
  const std::vector<int> &y_colind = input[1]->rowind();

  // loop over the rows of the resulting matrix)
  for(int i=0; i<z_rowind.size()-1; ++i){
    for(int el=z_rowind[i]; el<z_rowind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
      int j = z_col[el];
      int el1 = x_rowind[i];
      int el2 = y_colind[j];
      z_data[el] = 0;
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          z_data[el] |= x_data[el1++] | y_trans_data[el2++];
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
  }
}

} // namespace CasADi

