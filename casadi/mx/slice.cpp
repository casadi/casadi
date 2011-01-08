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

#include "slice.hpp"
#include <cassert>

using namespace std;

namespace CasADi{


Slice::Slice(const MX& x, Slicer i_, Slicer j_) : Reordering(x), i(i_), j(j_) {
  i.initialize(x.size1());
  j.initialize(x.size2());
  
  nrow_ = i.size();
  ncol_ = j.size();
}

Slice* Slice::clone() const{
  return new Slice(*this);
}

int Slice::k2k(int k) {
  return i((k/ncol_))*dep(0).size2() + j((k%ncol_));
}

void Slice::evaluate(int fsens_order, int asens_order){
 assert(fsens_order==0 || asens_order==0);
 
 


 
  if(fsens_order==0){
    int k=0;
    for (vector<int>::iterator iti=i.begin();iti!=i.end();++iti) {
     	for (vector<int>::iterator itj=j.begin();itj!=j.end();++itj) {
    	  val(0)[k++]=dep(0)->val(0)[(*iti)*dep(0).size2() + (*itj)];
      }
    }
  } else {
    int k=0;
    for (vector<int>::iterator iti=i.begin();iti!=i.end();++iti) {
     	for (vector<int>::iterator itj=j.begin();itj!=j.end();++itj) {
    	  val(1)[k++]=dep(0)->val(1)[(*iti)*dep(0).size2() + (*itj)];
      }
    }
  }
  
  if(asens_order>0){
    int k=0;
    for (vector<int>::iterator iti=i.begin();iti!=i.end();++iti) {
     	for (vector<int>::iterator itj=j.begin();itj!=j.end();++itj) {
    	  dep(0)->val(0)[(*iti)*dep(0).size2() + (*itj)]+=val(0)[k++];
      }
    }
  }
}

void Slice::print(std::ostream &stream) const{
  stream << "slice(" << dep(0) << ")";
}



} // namespace CasADi
