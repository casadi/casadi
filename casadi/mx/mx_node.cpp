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

MXNode::MXNode(const vector<MX>& dep) : dep_(dep){
  maxord_ = 0;
  nfdir_ = 1;
  nadir_ = 1;
}

MXNode::MXNode(const MX& dep){
  maxord_ = 0;
  nfdir_ = 1;
  nadir_ = 1;
  dep_.resize(1);
  dep_[0] = dep;
}

MXNode::MXNode(const MX& dep1, const MX& dep2){
  maxord_ = 0;
  nfdir_ = 1;
  nadir_ = 1;
  dep_.resize(2);
  dep_[0] = dep1;
  dep_[1] = dep2;
}

MXNode::MXNode(const MX& dep1, const MX& dep2, const MX& dep3){
  maxord_ = 0;
  nfdir_ = 1;
  nadir_ = 1;
  dep_.resize(3);
  dep_[0] = dep1;
  dep_[1] = dep2;
  dep_[2] = dep3;
}

MXNode::~MXNode(){
}

void MXNode::print(ostream &stream) const{
  stream << "<empty matrix expression>";
}

/*void MXNode::evaluateAdj(){
  cerr << "evaluateAdj not defined for class " << typeid(*this).name() << endl;
  throw "RuntimeElement::evaluateAdj";
}*/
  
const string& MXNode::getName() const{
  cerr << "getName not defined for class " << typeid(*this).name() << endl;
  throw "MXNode::getName()";
}

bool MXNode::isSymbolic() const{
  return false;
}

bool MXNode::isConstant() const{
  return false;
}

void MXNode::setOutput(const vector<double>& x, int ord){
  assert(x.size() == val(ord).size());
  copy(x.begin(),x.end(), val(ord).begin());
}

void MXNode::getOutput(vector<double>& x, int ord) const{
  assert(x.size() == val(ord).size());
  copy(val(ord).begin(),val(ord).end(),x.begin());
}

MX& MXNode::dep(int ind){
  return dep_.at(ind);
}

const MX& MXNode::dep(int ind) const{
  return dep_.at(ind);
}
  
const vector<double>& MXNode::val(int order, int dir) const{
  return val_.at(order).at(dir);
}

vector<double>& MXNode::val(int order, int dir){
  return val_.at(order).at(dir);
}

int MXNode::ndep() const{
  return dep_.size();
}

void MXNode::init(){
  val_.resize(maxord_ +1);
  for(int ord=0; ord<=maxord_; ++ord){
    val_[ord].resize(nfdir_+nadir_);
    for(int dir=0; dir<nfdir_+nadir_; ++dir){
      val_[ord][dir].resize(sz.nrow*sz.ncol);
    }
  }
}


} // namespace CasADi
