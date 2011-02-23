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

#include "generic_type_internal.hpp"
#include "stl_vector_tools.hpp"
#include "casadi_exception.hpp"

using namespace std;

namespace CasADi{

bool GenericTypeInternal::toBool() const{
  return bool(toInt());
}

int GenericTypeInternal::toInt() const{
  return toIntVector().at(0);
}

double GenericTypeInternal::toDouble() const{
  return toDoubleVector().at(0);
}

const vector<int>& GenericTypeInternal::toIntVector() const{
  if(n != i_vec.size()) throw CasadiException("GenericTypeInternal::toIntVector");
  return i_vec;
}

const vector<double>& GenericTypeInternal::toDoubleVector() const{
  if(n != d_vec.size()) throw CasadiException("GenericTypeInternal::toDoubleVector");
  return d_vec;
}

const string& GenericTypeInternal::toString() const{
  return str;
}

GenericTypeInternal::GenericTypeInternal(const vector<int>& i_vec_) : i_vec(i_vec_){
  stringstream ss;
  ss << i_vec;
  str = ss.str();
  is_string = true;
  n = i_vec.size();
  d_vec.resize(n);
  copy(i_vec.begin(),i_vec.end(),d_vec.begin());
}

GenericTypeInternal::GenericTypeInternal(const vector<double>& d_vec_) : d_vec(d_vec_){
  stringstream ss;
  ss << i_vec;
  str = ss.str();
  is_string = true;
  n = d_vec.size();
  i_vec.resize(n);
  copy(d_vec.begin(),d_vec.end(),i_vec.begin());
}

GenericTypeInternal::GenericTypeInternal(const string& s) : str(s){
  is_string = true;
  n = 1;
}


} // namespace CasADi

