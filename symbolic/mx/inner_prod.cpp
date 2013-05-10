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

#include "inner_prod.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../runtime/runtime.hpp"

using namespace std;

namespace CasADi{

  InnerProd::InnerProd(const MX& x, const MX& y){
    casadi_assert(x.sparsity()==y.sparsity());
    setDependencies(x,y);
    setSparsity(sp_dense(1,1));
  }
  
  void InnerProd::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "inner_prod(";
    } else if(part==1){
      stream << ",";
    } else {
      stream << ")";
    }
  }

  void InnerProd::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    if(!output_given){
      *output[0] = (*input[0])->getInnerProd(*input[1]);
    }

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for(int d=0; d<nfwd; ++d){
      *fwdSens[d][0] = (*input[0])->getInnerProd(*fwdSeed[d][1]) + (*fwdSeed[d][0])->getInnerProd(*input[1]);
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for(int d=0; d<nadj; ++d){
      *adjSens[d][0] += *adjSeed[d][0] * *input[1];
      *adjSens[d][1] += *adjSeed[d][0] * *input[0];
      *adjSeed[d][0] = MX();
    }
  }

  void InnerProd::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    evaluateGen<double,DMatrixPtrV,DMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  void InnerProd::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    evaluateGen<SX,SXMatrixPtrV,SXMatrixPtrVV>(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename T, typename MatV, typename MatVV>
  void InnerProd::evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens){
    // Get data
    T& res = output[0]->data().front();
    const vector<T> &arg0 = input[0]->data();
    const vector<T> &arg1 = input[1]->data();
    const int n = arg0.size();

    // Perform the inner product
    res = casadi_dot(n,getPtr(arg0),1,getPtr(arg1),1);

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for(int d=0; d<nfwd; ++d){
      T& fsens = fwdSens[d][0]->data().front();
      const vector<T> &fseed0 = fwdSeed[d][0]->data();
      const vector<T> &fseed1 = fwdSeed[d][1]->data();
      fsens = casadi_dot(n,getPtr(fseed0),1,getPtr(arg1),1) + casadi_dot(n,getPtr(arg0),1,getPtr(fseed1),1);
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for(int d=0; d<nadj; ++d){
      T& aseed = adjSeed[d][0]->data().front();
      vector<T> &asens0 = adjSens[d][0]->data();
      vector<T> &asens1 = adjSens[d][1]->data();
      casadi_axpy(n,aseed,getPtr(arg1),1,getPtr(asens0),1);
      casadi_axpy(n,aseed,getPtr(arg0),1,getPtr(asens1),1);
      aseed = 0;
    }
  }


} // namespace CasADi
