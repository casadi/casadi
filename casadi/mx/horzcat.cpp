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

#include "horzcat.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <iterator>

using namespace std;

namespace CasADi{

// Constructor
Horzcat::Horzcat(const vector<MX>& dep__) : MXNode(dep__){
  assert(!dep_.empty());
  int sz1=dep(0).size1();
  int sz2=0;
  cumulSize2.resize(dep_.size());
  int cnt=0;
  for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
    assert(sz1==it->size1());
    cumulSize2[cnt++] = sz2;
    sz2 += it->size2();
  }  
  nrow_ = sz1;
  ncol_ = sz2;
}

Horzcat* Horzcat::clone() const{
  return new Horzcat(dep_);
}

void Horzcat::print(ostream &stream) const{
  stream << "[";
  copy(dep_.begin(), dep_.end(), ostream_iterator<MX>(stream, ","));
  stream << "]";
}

void Horzcat::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || asens_order==0);
  
  if(fsens_order==0){
    int cnt=0;
    for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
      int offset = cumulSize2[cnt++];
      for (int k=0; k < it->size1(); k++) {
        copy(&(*it)->val(0)[k*it->size2()],&(*it)->val(0)[k*it->size2()]+it->size2(),&val(0)[offset+k*ncol_]);
      }
    }
  } else {
    int cnt=0;
    for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
      int offset = cumulSize2[cnt++];
      for (int k=0; k < it->size1(); k++) {
        copy(&(*it)->val(1)[k*it->size2()],&(*it)->val(1)[k*it->size2()]+it->size2(),&val(1)[offset+k*ncol_]);
      }
    }
  }
  
  if(asens_order>0){
    throw CasadiException("Horzcat::evaluate: adjoints evaluation not implemented");
  }
}

} // namespace CasADi
