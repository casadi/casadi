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

using namespace std;

namespace CasADi{

  InnerProd::InnerProd(const MX& x, const MX& y){
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
      casadi_error("not implemented");
      //*output[0] = det(*input[0]);
    }

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for(int d=0; d<nfwd; ++d){
      casadi_error("not implemented");
      //      *fwdSens[d][0] = trans(*fwdSeed[d][0]);
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for(int d=0; d<nadj; ++d){
      casadi_error("not implemented");
      //      *adjSens[d][0] +=  (*adjSeed[d][0]**output[0]) * trans(inv(*input[0]));
      //      *adjSeed[d][0] = MX();
    }
  }

} // namespace CasADi
