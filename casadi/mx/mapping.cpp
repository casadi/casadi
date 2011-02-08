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

#include "mapping.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

Mapping::Mapping(){
}

void Mapping::evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  
  for(int k=0; k<size(); ++k){
    output[k] = input[depind_[k]][nzind_[k]];
    
    for(int d=0; d<nfwd; ++d)
      fwdSens[d][k] = fwdSeed[depind_[k]][d][nzind_[k]];
    
    for(int d=0; d<nadj; ++d)
      adjSens[depind_[k]][d][nzind_[k]] += adjSeed[d][k];
  }
}

void Mapping::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "mapping(" << args << "," << depind_ << "," << nzind_ << ")" << endl;
}

void Mapping::addDepend(const MX& d, std::vector<int> nz, std::vector<int> i, std::vector<int> j){
  // Append the new dependency and save its index
  dep_.push_back(d);
  int k = dep_.size()-1;
  
  // New non-zero indices
  std::vector<int> nzind;
  nzind.reserve(nzind_.size());

  // New dependency index mapping
  std::vector<int> depind;
  depind.reserve(depind_.size());
 
  // Add the dependencies
  vector<int> el = sparsity_.getNZ(i,j);
  
  
  
  // Swap the vectors
  nzind_.swap(nzind);
  depind_.swap(depind);
}

} // namespace CasADi
