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

#include "reordering.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

Reordering::Reordering(){
}

void Reordering::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  
  for(int k=0; k<size(); ++k){
    output[k] = input[k2l(k)][k2k(k)];
    
    for(int d=0; d<nfwd; ++d)
      fwdSens[d][k] = fwdSeed[k2l(k)][d][k2k(k)];
    
    for(int d=0; d<nadj; ++d)
      adjSens[k2l(k)][d][k2k(k)] += adjSeed[d][k];
  }
}

int Reordering::k2l(int k) const{
  return argind_[k];
}

int Reordering::k2k(int k) const{
  return nzind_[k];
}

void Reordering::print(std::ostream &stream) const{
  stream << "reordering(" << dep(0) << "," << nzind_;
  if(ndep()>1) stream << "," << argind_;
  stream << ")";
}


} // namespace CasADi
