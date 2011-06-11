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

Slice::Slice(int start__, int stop__, int step__) : start(start__), stop(stop__), step(step__){
}

std::vector<int> Slice::getAll(int len) const{
  int start_ = start;
  int stop_  = stop;
  if (start_<0) start_+=len;
  if (stop_<=0) stop_+=len;
  return range(start_,stop_,step,len);
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
  } else if (type == SLICE) {
        return slice.getAll(len);
  }
}

  
} // namespace CasADi

