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
#include "../fx/fx_internal.hpp"

using namespace std;

namespace CasADi{

Multiplication::Multiplication(const MX& x, const MX& y_trans){
  casadi_assert_message(x.size2() == y_trans.size2(),"Multiplication::Multiplication: dimension mismatch. Attempting to multiply " << x.dimString() << " with " << y_trans.dimString());
  setDependencies(x,y_trans);

  // Create the sparsity pattern for the matrix-matrix product
  CRSSparsity spres = x->sparsity().patternProduct(y_trans.sparsity());

  // Save sparsity
  setSparsity(spres);
}

Multiplication* Multiplication::clone() const{
  return new Multiplication(*this);
}

void Multiplication::printPart(std::ostream &stream, int part) const{
  if(part==0){
    stream << "mul(";
  } else if(part==1){
    stream << ",trans(";
  } else {
    stream << "))";
  }
}

void Multiplication::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();

  fill(output[0]->begin(),output[0]->end(),0);
  DMatrix::mul_no_alloc_nt(*input[0],*input[1],*output[0]);

  // Forward sensitivities: dot(Z) = dot(X)*Y + X*dot(Y)
  for(int d=0; d<nfwd; ++d){
    fill(fwdSens[d][0]->begin(),fwdSens[d][0]->end(),0);
    DMatrix::mul_no_alloc_nt(*fwdSeed[d][0],*input[1],*fwdSens[d][0]);
    DMatrix::mul_no_alloc_nt(*input[0],*fwdSeed[d][1],*fwdSens[d][0]);
  }

  // Adjoint sensitivities
  for(int d=0; d<nadj; ++d){
    DMatrix::mul_no_alloc_nn(*adjSeed[d][0],*input[1],*adjSens[d][0]);
    DMatrix::mul_no_alloc_tn(*adjSeed[d][0],*input[0],*adjSens[d][1]);
    adjSeed[d][0]->setZero();
  }
}

void Multiplication::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  fill(output[0]->begin(),output[0]->end(),0);
  SXMatrix::mul_no_alloc_nt(*input[0],*input[1],*output[0]);
}

void Multiplication::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  if(!output_given)
    *output[0] = mul(*input[0],trans(*input[1]));

  // Forward sensitivities
  int nfwd = fwdSens.size();
  for(int d=0; d<nfwd; ++d){
    *fwdSens[d][0] = mul(*fwdSeed[d][0],trans(*input[1])) + mul(*input[0],trans(*fwdSeed[d][1]));
  }
  
  // Adjoint sensitivities
  int nadj = adjSeed.size();
  for(int d=0; d<nadj; ++d){
    *adjSens[d][0] += mul(*adjSeed[d][0],*input[1]);
    *adjSens[d][1] += mul(trans(*adjSeed[d][0]),*input[0]);
    *adjSeed[d][0] = MX();
  }
}

void Multiplication::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
  DMatrix::mul_sparsity(*input[0],*input[1],*output[0],fwd);
  if(!fwd) fill_n(get_bvec_t(output[0]->data()),output[0]->size(),bvec_t(0));
}

void Multiplication::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
  gen.addAuxiliary(CodeGenerator::AUX_MM_NT_SPARSE);
  gen.addAuxiliary(CodeGenerator::AUX_FILL);

  // Clear the result
  stream << "  casadi_fill(" << sparsity().size() << ",0.0," << res.front() << ",1);" << endl;

  // Perform sparse matrix multiplication
  stream << "  casadi_mm_nt_sparse(";
  for(int i=0; i<2; ++i){
    stream << arg.at(i) << ",s" << gen.getSparsity(dep(i).sparsity()) << ",";
  }
  stream << res.front() << ",s" << gen.getSparsity(sparsity()) << ");" << endl;
}


} // namespace CasADi

