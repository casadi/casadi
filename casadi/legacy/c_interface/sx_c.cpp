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

#include "sx_c.hpp"
#include "stl_string_c.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;
namespace CasADi{
  
SX& sx_ref(sx_ptr ptr){
  if(ptr==0) throw CasadiException("get_sx failed: null pointer");
  SX *r = (SX*)(ptr);
  return *r;
}

std::vector<SX>& sx_vec(sx_vec_ptr ptr){
  if(ptr==0) throw CasadiException("sx_vec failed: null pointer");
  vector<SX> *r = (vector<SX>*)(ptr);
  return *r;
}

extern "C"
sx_ptr casadi_sx_new(){
  c_constructor_wrapper("casadi_sx_new",
    return new SX();
  )
}

extern "C"
int casadi_sx_symbol(sx_ptr ptr, const char* name){
  c_fcn_wrapper("casadi_sx_symbol",
    sx_ref(ptr) = SX(name);
  )
}

extern "C"
int casadi_sx_constant(sx_ptr ptr, double value){
  c_fcn_wrapper("casadi_sx_constant",
    sx_ref(ptr) = SX(value);
  )
}

extern "C"
int casadi_sx_delete(sx_ptr ptr){
  c_fcn_wrapper("casadi_sx_delete",
    delete &sx_ref(ptr);
  )
}

extern "C"
int casadi_sx_print(sx_ptr ptr, string_ptr str){
  c_fcn_wrapper("casadi_sx_print_string",
    // Print to a temporary variable
    stringstream ss;
    ss << sx_ref(ptr);
    get_stl_string(str) = ss.str();
  )
}

extern "C"
int casadi_sx_binary(sx_ptr r, int op, sx_ptr x, sx_ptr y){
  c_fcn_wrapper("casadi_sx_binary",
    sx_ref(r) = SX::binary(op,sx_ref(x),sx_ref(y));
  )
}

extern "C"
int casadi_sx_unary(sx_ptr r, int op, sx_ptr x){
  c_fcn_wrapper("casadi_sx_binary",
    sx_ref(r) = SX::unary(op,sx_ref(x));
  )
}

extern "C"
sx_vec_ptr casadi_sx_vec_new(void){
  c_constructor_wrapper("casadi_sx_vec_new",
    return new vector<SX>();
  )
}

extern "C"
int casadi_sx_vec_delete(sx_vec_ptr ptr){
  c_fcn_wrapper("casadi_sx_delete",
    delete &sx_vec(ptr);
  )
}

extern "C"
int casadi_sx_vec_push_back(sx_vec_ptr v, sx_ptr ptr){
  c_fcn_wrapper("casadi_sx_vec_push_back",
    sx_vec(v).push_back(sx_ref(ptr));
  )
}

extern "C"
int casadi_sx_vec_size(sx_vec_ptr v){
  try{
    return sx_vec(v).size();
  } catch(...){
    return -1;
  }
}

extern "C"
int casadi_sx_vec_print(sx_vec_ptr ptr, string_ptr str){
  c_fcn_wrapper("casadi_sx_vec_print",
    // Print to a temporary variable
    stringstream ss;
    ss << sx_vec(ptr);
    get_stl_string(str) = ss.str();
  )
}


} // namespace CasADi
