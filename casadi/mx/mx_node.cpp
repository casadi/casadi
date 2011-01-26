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

#include "mx_node.hpp"
#include <cassert>
#include <typeinfo> 

using namespace std;

namespace CasADi{

MXNode::MXNode(){
}

MXNode::~MXNode(){
}

const string& MXNode::getName() const{
  throw CasadiException(string("MXNode::getName() not defined for class ") + typeid(*this).name());
}

bool MXNode::isSymbolic() const{
  return false;
}

bool MXNode::isConstant() const{
  return false;
}

const MX& MXNode::dep(int ind) const{
  return dep_.at(ind);
}
  
int MXNode::ndep() const{
  return dep_.size();
}

void MXNode::setSparsity(const CRSSparsity& sparsity){
  sparsity_ = sparsity;
}

void MXNode::setDependencies(const MX& dep){
  dep_.resize(1);
  dep_[0] = dep;
}
    
void MXNode::setDependencies(const MX& dep1, const MX& dep2){
  dep_.resize(2);
  dep_[0] = dep1;
  dep_[1] = dep2;
}
    
void MXNode::setDependencies(const MX& dep1, const MX& dep2, const MX& dep3){
  dep_.resize(3);
  dep_[0] = dep1;
  dep_[1] = dep2;
  dep_[2] = dep3;
}

void MXNode::setDependencies(const std::vector<MX>& dep){
  dep_ = dep;
}

int MXNode::numel() const{
  return sparsity_.numel();
}

int MXNode::size() const{
  return sparsity_.size();
}

int MXNode::size1() const{
  return sparsity_.size1();
}

int MXNode::size2() const{
  return sparsity_.size2();
}

const CRSSparsity& MXNode::sparsity() const{
  return sparsity_;
}

void MXNode::print(std::ostream &stream) const{
  vector<string> args(ndep());
  for(int i=0; i<ndep(); ++i){
    stringstream ss;
    dep(i)->print(ss);
    args[i] = ss.str();
  }
  print(stream,args);
}

std::vector<MX> MXNode::partial() const{
  throw CasadiException(string("MXNode::partial() no partial derivatives defined for class ") + typeid(*this).name());
}



} // namespace CasADi
