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
  setSize(i.size(),j.size());
}

Slice* Slice::clone() const{
  return new Slice(*this);
}

void Slice::init(){
  // Base class initializations
  Reordering::init();
  
  // Mapping: NOTE: does not handle sparsity correctly, should get the reordering from the Matrix<> class
  for(int k=0; k<size1()*size2(); ++k)
    nzind_[k] = i((k/size2()))*dep(0).size2() + j((k%size2()));
}

void Slice::evaluate(int fsens_order, int asens_order){
 assert(fsens_order==0 || asens_order==0);
 
  if(fsens_order==0){
    int k=0;
    for (vector<int>::iterator iti=i.begin();iti!=i.end();++iti) {
     	for (vector<int>::iterator itj=j.begin();itj!=j.end();++itj) {
    	  output()[k++]=input(0)[(*iti)*dep(0).size2() + (*itj)];
      }
    }
  } else {
    int k=0;
    for (vector<int>::iterator iti=i.begin();iti!=i.end();++iti) {
     	for (vector<int>::iterator itj=j.begin();itj!=j.end();++itj) {
    	  fwdSens()[k++]=fwdSeed(0)[(*iti)*dep(0).size2() + (*itj)];
      }
    }
  }
  
  if(asens_order>0){
    int k=0;
    for (vector<int>::iterator iti=i.begin();iti!=i.end();++iti) {
     	for (vector<int>::iterator itj=j.begin();itj!=j.end();++itj) {
          input(0)[(*iti)*dep(0).size2() + (*itj)]+=output()[k++]; // NOTE TO JORIS: This must be wrong
      }
    }
  }
}

void Slice::print(std::ostream &stream) const{
  stream << "slice(" << dep(0) << ")";
}



} // namespace CasADi
