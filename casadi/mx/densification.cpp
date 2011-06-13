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
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

Densification::Densification(const MX& x){
  setDependencies(x);
  
  // Calculate the mapping of the nonzeros
  CRSSparsity sp = x->sparsity().makeDense(mapping_);
  setSparsity(sp);
}

Densification* Densification::clone() const{
  return new Densification(*this);
}

void Densification::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "dense(" << args.at(0) << ")";
}

void Densification::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  for(int el=0; el<mapping_.size(); ++el)
    output[mapping_[el]] = input[0][el];
  
  // Propagate forward seeds
  for(int d=0; d<nfwd; ++d){
    for(int el=0; el<mapping_.size(); ++el){
      fwdSens[d][mapping_[el]] = fwdSeed[0][d][el];
    }
  }
      
  // Propagate adjoint seeds
  for(int d=0; d<nadj; ++d){
    for(int el=0; el<mapping_.size(); ++el){
      adjSens[0][d][el] = adjSeed[d][mapping_[el]];
    }
  }
}

MX Densification::adFwd(const std::vector<MX>& jx){
  return jx.at(0);
}



} // namespace CasADi

