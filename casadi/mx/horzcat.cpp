/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
  int sz1=dep_[0].size1();
  int sz2=0;
  for(vector<MX>::const_iterator it=dep_.begin(); it!=dep_.end(); ++it){
    assert(sz1==it->size1());
    sz2 += it->size2();
  }  
  sz.nrow = sz1;
  sz.ncol = sz2;
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
  assert(0);
}

} // namespace CasADi
