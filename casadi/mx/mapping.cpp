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
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace CasADi{

Mapping::Mapping(const CRSSparsity& sp){
  setSparsity(sp);
  nzind_.resize(sp.size(),-1);
  depind_.resize(sp.size(),-1);
}

Mapping* Mapping::clone() const{
  return new Mapping(*this);
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

bool Mapping::isReady() const{
  casadi_assert(depind_.size()==size());
  casadi_assert(nzind_.size()==size());
  for(int k=0; k<size(); ++k){
    if(nzind_[k]<0 || depind_[k]<0)
      return false;
  }
  return true;
}
    

void Mapping::print(std::ostream &stream, const std::vector<std::string>& args) const{
  casadi_assert(isReady());
  
  if(numel()==1 && size()==1 && ndep()==1){
    stream << args[0];
    if(dep(0).numel()>1) stream << "[" << nzind_.at(0) << "]";
  }
  else{
    stream << "mapping(" << size1() << "-by-" << size2() << " matrix, nonzeros: [";
    for(int i=0; i<nzind_.size(); ++i){
      if(i!=0) stream << ",";
      stream << args[depind_[i]];
      if(dep(depind_[i]).numel()>1) stream << "[" << nzind_[i] << "]";
    }
    stream << "])";
  }
}

void Mapping::addDependency(const MX& d, const std::vector<int>& nz_d){
  addDependency(d,nz_d,range(nz_d.size()));
}

void Mapping::addDependency(const MX& d, const std::vector<int>& nz_d, const std::vector<int>& nz){
  casadi_assert(nz_d.size()==nz.size());
  
  // Quick return if no elements
  if(nz_d.empty()) return;
  
  if(true && d->isMapping()){
    // Eliminate if a mapping node
    const Mapping* dnode = static_cast<const Mapping*>(d.get());
    vector<MX> d2 = dnode->dep_;
    vector<vector<int> > nz_d2(d2.size());
    vector<vector<int> > nz2(d2.size());
    for(int i=0; i<nz.size(); ++i){
      int depind_i = dnode->depind_.at(nz_d[i]);
      nz_d2[depind_i].push_back(dnode->nzind_.at(nz_d[i]));
      nz2[depind_i].push_back(nz[i]);
    }
    
    // Call the function recursively
    for(int i=0; i<d2.size(); ++i){
      addDependency(d2[i],nz_d2[i],nz2[i]);
    }
  } else {
    // Add the node if it is not already a dependency
    std::map<const MXNode*, int>::const_iterator it = depmap_.find(static_cast<const MXNode*>(d.get()));
    int depind;
    if(it==depmap_.end()){
      depind = MXNode::addDependency(d);
      depmap_[static_cast<const MXNode*>(d.get())] = depind;
    } else {
      depind = it->second;
    }
    
    // Save the mapping
    addDependency(depind,nz_d,nz);
  }
}

void Mapping::addDependency(int depind, const std::vector<int>& nz_d, const std::vector<int>& nz){
  casadi_assert(nz_d.size()==nz.size());
  for(int k=0; k<nz.size(); ++k){
    nzind_[nz[k]] = nz_d[k];
    depind_[nz[k]] = depind;
  }
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
