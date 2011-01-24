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

#include "reshape.hpp"
#include <algorithm>

using namespace std;

namespace CasADi{

Reshape::Reshape(const MX& x, int n, int m){
  if (n*m != x.numel()) {
    throw CasadiException("MX::reshape: size must be same before and after reshaping");
  }
  setDependencies(x);
  setSparsity(CRSSparsity(n,m,true));
}

Reshape* Reshape::clone() const{
  return new Reshape(*this);
}

void Reshape::print(std::ostream &stream) const{
  stream << "reshape(" << dep(0) << ",[" << size1() << "," << size2() << "])";
}

void Reshape::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  for(int i=0; i<size(); ++i){
    // Function
    output[i] = input[0][i];
    
    // Forward seeds
    for(int d=0; d<nfwd; ++d)
      fwdSens[d][i] = fwdSeed[0][d][i];
    
    // Adjoint seeds
    for(int d=0; d<nadj; ++d)
      adjSens[0][d][i] += adjSeed[d][i];
  }
}

} // namespace CasADi
