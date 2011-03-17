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

#include "mx_tools.hpp"
#include "mapping.hpp"
#include "norm.hpp"
#include "multiplication.hpp"
#include "if_else_node.hpp"
#include "../fx/mx_function.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

MX vertcat(const vector<MX>& comp){
  // Remove nulls
  vector<MX> c;
  c.reserve(comp.size());
  for(vector<MX>::const_iterator it=comp.begin(); it!=comp.end(); ++it)
    if(!it->isNull())
      c.push_back(*it);
  
  if(c.empty()){
    return MX();
  } else if(c.size()==1){
    return c[0];
  } else {
    // Construct the sparsity pattern
    CRSSparsity sp = c[0].sparsity();
    for(int i=1; i<c.size(); ++i){
      sp.append(c[i].sparsity());
    }

    // Create a mapping matrix with the corresponding sparsity
    MX ret;
    ret.assignNode(new Mapping(sp));
    
    // Map the dependencies
    int offset=0;
    for(int i=0; i<c.size(); ++i){
      int nz = c[i].size();
      ret->addDependency(c[i],range(nz),range(offset,offset+nz));
      offset += nz;
    }
    
    return ret;
  }
}

MX horzcat(const vector<MX>& comp){
  vector<MX> v(comp.size());
  for(int i=0; i<v.size(); ++i)
    v[i] = trans(comp[i]);
  return trans(vertcat(v));
}

MX vertcat(const MX& a, const MX& b){
  vector<MX> ab;
  ab.push_back(a);
  ab.push_back(b);
  return vertcat(ab);
}

MX horzcat(const MX& a, const MX& b){
  vector<MX> ab;
  ab.push_back(a);
  ab.push_back(b);
  return horzcat(ab);
}

MX norm_2(const MX &x){
  MX ret;
  ret.assignNode(new Norm2(x));
  return ret;
}

MX norm_1(const MX &x){
  MX ret;
  ret.assignNode(new Norm1(x));
  return ret;
}

MX norm_inf(const MX &x){
  MX ret;
  ret.assignNode(new NormInf(x));
  return ret;
}

MX prod(const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new Multiplication(x,y));
  return ret;
}

MX inner_prod(const MX &x, const MX &y){
  return prod(trans(x),y);
}

MX outer_prod(const MX &x, const MX &y){
  return prod(x,trans(y));
}

MX trans(const MX &x){
  // Quick return if null or scalar
  if(x.isNull() || x.numel()==1)
    return x;
  
  // Get the tranposed matrix and the corresponding mapping
  vector<int> nzind;
  CRSSparsity sp = x->sparsity().transpose(nzind);

  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(x,nzind);
  return ret;
}

MX reshape(const MX &x, const std::vector<int> sz){
  if(sz.size() != 2)
    throw CasadiException("MX::reshape: not two dimensions");
  return reshape(x,sz[0],sz[1]);
}

MX reshape(const MX &x, int n, int m){
  // only works if dense
  casadi_assert_message(x.size()==x.numel(),"Sparse reshape not implemented");
  return reshape(x,CRSSparsity(n,m,true));
}

MX reshape(const MX &x, const CRSSparsity& sp){
  // quick return if already the right shape
  if(sp==x.sparsity())
    return x;
  
  // make sure that the number of zeros agree
  casadi_assert(x.size()==sp.size());
  
  // Create a mapping
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(x,range(x.size()));
  return ret;
}

MX vec(const MX &x) {
  // only works if dense
  return reshape(x,x.numel(),1);
}

MX if_else_zero(const MX &cond, const MX &if_true){
  MX ret;
  ret.assignNode(new IfNode(cond,if_true));
  return ret;
}


MX if_else(const MX &cond, const MX &if_true, const MX &if_false){
  return if_else_zero(cond,if_true) + if_else_zero(1-cond,if_false);
}

MX unite(const MX& A, const MX& B){
  // Join the sparsity patterns
  std::vector<int> mapping;
  CRSSparsity sp = A.sparsity().patternUnion(B.sparsity(),mapping);
  
  // Split up the mapping
  std::vector<int> nzA,nzB;
  
  // Copy sparsity
  for(int k=0; k<mapping.size(); ++k){
    if(mapping[k]<0){
      nzA.push_back(k);
    } else if(mapping[k]>0){
      nzB.push_back(k);
    } else {
      throw CasadiException("Pattern intersection not empty");
    }
  }
  
  // Create mapping
  MX ret;
  ret.assignNode(new Mapping(sp));
  ret->addDependency(A,range(nzA.size()),nzA);
  ret->addDependency(B,range(nzB.size()),nzB);
  return ret;
}



} // namespace CasADi

