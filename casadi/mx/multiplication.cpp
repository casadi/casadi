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

void Multiplication::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
  fill(output.begin(),output.end(),0);
  SXMatrix::prod_no_alloc(*input[0],*input[1],output);
}

MX Multiplication::adFwd(const std::vector<MX>& jx){
  
  // Get the number of derivative directions
  int nd = jx[0].size2();
  
  // The columns of the jacobian
  vector<MX> jac_cols;
  
  // Loop over the derivative directions
  for(int d=0; d<nd; ++d){
    
    // Get a forward seed matrices for direction d
    MX seed0 = reshape(jx[0](range(jx[0].size1()),d),dep(0).size1(),dep(0).size2());
    MX seed1 = reshape(jx[1](range(jx[1].size1()),d),dep(1).size1(),dep(1).size2());

    // Calculate the forward sensitivity
    MX sens = prod(seed0,trans(dep(1))) + prod(dep(0),trans(seed1));
   
    // Save the column of the Jacobian
    jac_cols.push_back(vec(sens));
  }
  
  // Assemble all the directions and return
  return horzcat(jac_cols);
}



} // namespace CasADi

