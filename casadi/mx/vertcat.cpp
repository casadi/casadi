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
#include <iterator>
#include <algorithm>

using namespace std;

namespace CasADi{

// Constructor
Vertcat::Vertcat(const vector<MX>& d){
  setDependencies(d);
  if(d.empty())
    throw CasadiException("Vertcat: empty concatenation not allowed");
  int sz1=0;
  int sz2=dep(0).size2();
  for(int i=0; i<ndep(); ++i){
    sz1 += dep(i).size1();
    if(sz2!=dep(i).size2())
      throw CasadiException("Vertcat: dimension mismatch");
  }
  setSize(sz1,sz2);
}

Vertcat* Vertcat::clone() const{
  return new Vertcat(*this);
}

void Vertcat::print(ostream &stream) const{
  stream << "[";
  copy(dep_.begin(), dep_.end(), ostream_iterator<MX>(stream, ";"));
  stream << "]";
}

void Vertcat::evaluate(int fsens_order, int asens_order){
  if(fsens_order==0){
    vector<double>::iterator out = output().begin();
    for(int i=0; i<ndep(); ++i){
      copy(input(i).begin(),input(i).end(),out);
      out += input(i).size();
    }
  } else {
    vector<double>::iterator fsens = fwdSens().begin();
    for(int i=0; i<ndep(); ++i){
      copy(fwdSeed(i).begin(),fwdSeed(i).end(),fsens);
      fsens += fwdSeed(i).size();
    }
  }
  
  if(asens_order>0){
    int offset = 0;
    for(int i=0; i<ndep(); ++i){
      transform(adjSeed().begin()+offset,adjSeed().begin()+offset+adjSens(i).size(),
                adjSens(i).begin(),
                adjSens(i).begin(),
                plus<double>());
      offset += adjSens(i).size();
    }
  }
}


} // namespace CasADi
