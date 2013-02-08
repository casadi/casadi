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
#include "matrix_tools.hpp"

namespace CasADi{

Slice::Slice() : start_(0), stop_(std::numeric_limits<int>::max()), step_(1){ 
}
      
Slice::Slice(int i) : start_(i), stop_(i+1), step_(1){
}

Slice::Slice(int start, int stop, int step) : start_(start), stop_(stop), step_(step){
}

std::vector<int> Slice::getAll(int len) const{
  int start = start_;
  int stop  = stop_;
  if (start<0) start+=len;
  if (stop<0) stop+=len;
  if (stop==std::numeric_limits<int>::max()) stop = len;

  casadi_assert_message(stop<=len,"Slice (start=" << start << ", stop=" << stop << ", step=" << step_ << ") out of bounds with supplied length of " << len);
  casadi_assert_message(start>=0, "Slice (start=" << start << ", stop=" << stop << ", step=" << step_ << ") out of bounds with start<0.");
  if (stop>=start && step_<0 || stop<=start && step_>0) return std::vector<int>();

  return range(start,stop,step_,len);
}

IndexSet::IndexSet(int i) : v_(1,i){
}
  
IndexSet::IndexSet(const std::vector<int>& v) : v_(v){
}
      
const std::vector<int>& IndexSet::getAll(int len) const{
  return v_;
}


IndexList::IndexList() : type(NILL) {}
IndexList::IndexList(int i_) : i(i_), type(INT) {}
IndexList::IndexList(const std::vector<int> &i_) : iv(i_) ,type(IVECTOR) {}
IndexList::IndexList(const Slice &s) : slice(s), type(SLICE) {}

std::vector<int> IndexList::getAll(int len) const {
  if (type == INT)  {
    if (i<0) return std::vector<int>(1,i+len);
    return std::vector<int>(1,i);
  } else if (type == IVECTOR) {
    return iv;
  } else {
    casadi_assert(type == SLICE);
    return slice.getAll(len);
  }
}

  
} // namespace CasADi

