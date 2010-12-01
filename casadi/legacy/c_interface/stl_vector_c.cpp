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

#include "stl_vector_c.hpp"

using namespace std;

std::vector<double>& get_stl_vector(vector_ptr v){
  if(v==0) throw "get_stl_vector failed: null pointer";
  vector<double> *r = (vector<double>*)(v);
  return *r;
}

extern "C"
vector_ptr casadi_vector_new(void){
  try{
    return new vector<double>();
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_vector_delete(vector_ptr v){
  try{
    vector<double> *ptr = &get_stl_vector(v);
    delete ptr;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_vector_size(vector_ptr v){
  try{
    return get_stl_vector(v).size();
  } catch(...){
    return -1;
  }
}

extern "C"
const double* casadi_vector_get_ptr(vector_ptr v){
  try{
    return &get_stl_vector(v)[0];
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_vector_resize(vector_ptr v, int len){
  try{
    get_stl_vector(v).resize(len);
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_vector_set(vector_ptr v, const double* val){
  try{
    vector<double>& vec = get_stl_vector(v);
    copy(val,val+vec.size(),vec.begin());
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_vector_get(vector_ptr v, double* val){
  try{
    const vector<double>& vec = get_stl_vector(v);
    copy(vec.begin(),vec.end(),val);
  } catch(...){
    return 1;
  }  
}

