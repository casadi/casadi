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

#include "vertcat.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <iterator>

using namespace std;

namespace CasADi{

// Constructor
Vertcat::Vertcat(const vector<MX>& dep__) : MXNode(dep__){
  assert(!dep_.empty());
  int sz1=0;
  int sz2=dep(0).size2();
  for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
    sz1 += it->size1();
    assert(sz2==it->size2());
  }
  setSize(sz1,sz2);
}

Vertcat* Vertcat::clone() const{
  return new Vertcat(dep_);
}

void Vertcat::print(ostream &stream) const{
  stream << "[";
  copy(dep_.begin(), dep_.end(), ostream_iterator<MX>(stream, ";"));
  stream << "]";
}

void Vertcat::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
    int i = 0;
    for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
      copy((*it)->output().begin(),(*it)->output().end(),&output()[i]);
      i += it->numel();
    }
  } else {
    int i = 0;
    for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
      copy((*it)->fwdSens().begin(),(*it)->fwdSens().end(),&fwdSens()[i]);
      i += it->numel();
    }
  }
  
  if(asens_order>0){
    int i = 0;
    for(vector<MX>::iterator it=dep_.begin(); it!=dep_.end(); ++it){
      copy(&adjSeed()[i],&adjSeed()[i] + it->numel(), (*it)->adjSeed().begin());
      i += it->numel();
    }
  }
}

// void Vertcat::setOutput(const vector<double>& x, int ord){
//   int i = 0;
//   for(vector<MX>::iterator it=dep_.begin(); it!=dep_.end(); ++it){
//     (*it)->setOutput(&x[i],ord);
//     i += it->size();
//   }
// }
// 
// void Vertcat::getOutput(vector<double>& x, int ord) const{
//   int i = 0;
//   for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
//     (*it)->getOutput(&x[i],ord);
//     i += it->size();
//   }
// }



} // namespace CasADi
