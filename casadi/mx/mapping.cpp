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

Mapping::Mapping(const CRSSparsity& sp) : nzmap_(sp,-1){
  setSparsity(sp);
  depind_.resize(sp.size(),-1);
}

Mapping* Mapping::clone() const{
  return new Mapping(*this);
}

void Mapping::evaluate(const std::vector<DMatrix*>& input, DMatrix& output, const VVDptr& fwdSeed, std::vector<DMatrix*>& fwdSens, const std::vector<DMatrix*>& adjSeed, VVDptr& adjSens, int nfwd, int nadj){
  const std::vector<int>& nzind_ = nzmap_.data();
  vector<double> &outputd = output.data();
  
  for(int k=0; k<size(); ++k){
    outputd[k] = input[depind_[k]]->data()[nzind_[k]];
    
    for(int d=0; d<nfwd; ++d)
      fwdSens[d]->data()[k] = fwdSeed[depind_[k]][d][nzind_[k]];
    
    for(int d=0; d<nadj; ++d)
      adjSens[depind_[k]][d][nzind_[k]] += adjSeed[d]->data()[k];
  }
}

bool Mapping::isReady() const{
  const std::vector<int>& nzind_ = nzmap_.data();
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
  const std::vector<int>& nzind_ = nzmap_.data();
  
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
  casadi_assert(!d.isNull());
  const std::vector<int>& nzind_ = nzmap_.data();
  
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
      nz_d2[depind_i].push_back(dnode->nzmap_.at(nz_d[i]));
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
  std::vector<int>& nzind_ = nzmap_.data();
  for(int k=0; k<nz.size(); ++k){
    nzind_[nz[k]] = nz_d[k];
    depind_[nz[k]] = depind;
  }
}

MX Mapping::adFwd(const std::vector<MX>& jx){
  casadi_assert(isReady());
  casadi_assert(size()==numel());
  const std::vector<int> &nzind_ = nzmap_.data();

  // Number of columns
  int ncol = jx.front().size2();
    
  // Sparsity
  const CRSSparsity &sp = sparsity();

  // Nonzero elements of the new matrix
  vector<int> i_ret, j_ret, el_ret, dp_ret;

  // Loop over rows of the matrix
  for(int i=0; i<size1(); ++i){
    
    // Loop over nonzero entries
    for(int el=sp.rowind(i); el<sp.rowind(i+1); ++el){
      
      // Get the column
      //int j = sp.col(el);
      
      // Get the dependency and nonzero to which the nonzero is mapped
      int dp = depind_[el];
      int nz = nzind_[el];
      
      // Get the sparsity of the seed matrix
      const CRSSparsity &sp_mat = jx[dp].sparsity();
      
      // Loop over the nonzeros of the row of the seed matrix
      int i_mat = nz;
      for(int el_mat=sp_mat.rowind(i_mat); el_mat<sp_mat.rowind(i_mat+1); ++el_mat){
        // Get the column
        int j_mat = sp_mat.col(el_mat);
        
        // Map the Jacobian nonzero
        i_ret.push_back(el);
        j_ret.push_back(j_mat);
        el_ret.push_back(el_mat);
        dp_ret.push_back(dp);
      }
    }
  }
  
  // Row offsets for the return matrix
  vector<int> rowind(1,0);
  rowind.reserve(size()+1);
  int i=0;
  for(int k=0; k<i_ret.size(); ++k){
    casadi_assert(i_ret[k]>=i);
    for(; i<i_ret[k]; ++i){
      rowind.push_back(k);
    }
  }
  rowind.resize(size()+1,i_ret.size());
  
  // Sparsity of the return matrix
  CRSSparsity sp_ret(size(),ncol,j_ret,rowind);
  
  // Return matrix
  MX ret = MX::create(new Mapping(sp_ret));
  
  // Add the dependencies
  for(int dp=0; dp<jx.size(); ++dp){
    
    // Get the local nonzeros
    vector<int> nz, nzd;
    for(int k=0; k<el_ret.size(); ++k){
      
      // If dependency matches
      if(dp_ret[k]==dp){
        nz.push_back(k);
        nzd.push_back(el_ret[k]);
      }
    }
    
    // Save to return matrix
    ret->addDependency(jx[dp],nzd,nz);
  }
    
  // Return the mapping
  return ret;
}

void Mapping::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
  const std::vector<int> &nzind_ = nzmap_.data();
  for(int k=0; k<size(); ++k){
    output[k] = (*input[depind_[k]])[nzind_[k]];
  }
}

// MX Mapping::eval(const std::vector<MX>& x){
//   
// }

} // namespace CasADi
