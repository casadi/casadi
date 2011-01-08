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
#include "vertcat.hpp"
#include "horzcat.hpp"
#include "transpose.hpp"
#include "flatten.hpp"
#include "reshape.hpp"
#include "norm.hpp"
#include "multiplication.hpp"


namespace CasADi{

MX vertcat(const vector<MX>& comp){
  MX ret;
  ret.assignNode(new Vertcat(comp));
  return ret;
}

MX horzcat(const vector<MX>& comp){
  MX ret;
  ret.assignNode(new Horzcat(comp));
  return ret;
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
  // Check if the node is already a transpose
  const Transpose* t = dynamic_cast<const Transpose*>(x.get());

  if(t) // already a transpose
    return t->dep(0);
  else{
    MX ret;
    ret.assignNode(new Transpose(x));
    return ret;
  }
}

MX reshape(const MX &x, const std::vector<int> sz){
  if(sz.size() != 2)
    throw CasadiException("MX::reshape: not two dimensions");
  return reshape(x,sz[0],sz[1]);
}

MX reshape(const MX &x, int n, int m){
  if (n*m!=x.numel()) {
    throw CasadiException("MX::reshape: size must be same before and after reshaping");
  }
  if (n==x.size1() && m==x.size2()) {
    // allready correct shape
    return x.get()->dep(0);
  } else {
    MX ret;
    ret.assignNode(new Reshape(x,n,m));
    return ret;
  }
        
}

MX flatten(const MX &x) {
  if (x.size2()==1) {
          // Allready flattened
          return x.get()->dep(0);
  } else {
    MX ret;
    ret.assignNode(new Flatten(x));
    return ret;
  }
}


  
} // namespace CasADi

