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

#include "determinant.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"

using namespace std;

namespace CasADi{

  Determinant::Determinant(const MX& x){
    setDependencies(x);
    setSparsity(sp_dense(1,1));
  }
  
  void Determinant::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "det(";
    } else {
      stream << ")";
    }
  }

  void Determinant::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    int nfwd = fwdSens.size();    
    int nadj = adjSeed.size();

    // Non-differentiated output
    const MX& X = *input[0];
    MX& det_X = *output[0];
    if(!output_given){
      det_X = det(X);
    }
    
    // Quick return
    if(nfwd==0 && nadj==0) return;
    
    // Create only once
    MX trans_inv_X = trans(inv(X));

    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      *fwdSens[d][0] = det_X * inner_prod(trans_inv_X, *fwdSeed[d][0]);
    }

    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      *adjSens[d][0] +=  (*adjSeed[d][0]*det_X) * trans_inv_X;
      *adjSeed[d][0] = MX();
    }
  }

} // namespace CasADi
