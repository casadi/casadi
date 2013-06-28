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

#include "inverse.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"

using namespace std;

namespace CasADi{

  Inverse::Inverse(const MX& x){
    setDependencies(x);
    setSparsity(sp_dense(x.size1(),x.size2()));
  }
  
  void Inverse::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "inv(";
    } else {
      stream << ")";
    }
  }

  void Inverse::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    const MX& X = *input[0];
    MX& inv_X = *output[0];
    if(!output_given){
      inv_X = inv(X);
    }

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for(int d=0; d<nfwd; ++d){
      *fwdSens[d][0] = -mul(mul(inv_X,*fwdSeed[d][0]),inv_X);
    }
  
    // Adjoint sensitivities
    int nadj = adjSeed.size();
    if(nadj>0){
      MX trans_inv_X = trans(inv_X);
      for(int d=0; d<nadj; ++d){
        *adjSens[d][0] -= mul(mul(trans_inv_X,*adjSeed[d][0]),trans_inv_X);
        *adjSeed[d][0] = MX();
      }
    }
  }

} // namespace CasADi
