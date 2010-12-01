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

#include "sx_matrix_c.hpp"
#include "stl_string_c.hpp"
#include "../fx/sx_function.hpp"
#include "../stl_vector_tools.hpp"
#include "../expression_tools.hpp"

#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;
namespace CasADi{

// Get a reference
SXMatrix& get_sx_matrix(sx_matrix_ref ref){
  if(ref==0) throw CasadiException("get_sx_matrix failed: null pointer");
  SXMatrix *r = (SXMatrix*)(ref);
  return *r;
}

// Get a reference to a vector
vector<SXMatrix>& get_sx_matrix_vec(sx_matrix_vec v){
  if(v==0) throw CasadiException("get_sx_matrix_vec failed: null pointer");
  vector<SXMatrix> *r = (vector<SXMatrix>*)(v);
  return *r;
}

extern "C"
sx_matrix_ref casadi_sx_matrix_new(){
  return new SXMatrix();
}

extern "C"
int casadi_sx_matrix_symbol(sx_matrix_ref ref, const char* name, int nrow, int ncol){
  try{
    get_sx_matrix(ref) = SXMatrix(name,nrow,ncol);
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_constant(sx_matrix_ref ref, const double* data, int nrow, int ncol, char order){
  try{
    vector<double> temp(data,data+nrow*ncol);
    get_sx_matrix(ref) = SXMatrix(temp,nrow,ncol);
    return 0;
  } catch(...){
    return 1;
  }  
}

extern "C"
int casadi_sx_matrix_delete(sx_matrix_ref ref){
  SXMatrix* ptr = &get_sx_matrix(ref);
  delete ptr;
  return 0;
}

extern "C"
int casadi_sx_matrix_print(sx_matrix_ref ref, string_ptr str){
  c_fcn_wrapper("casadi_sx_matrix_print",
    stringstream ss;
    ss << get_sx_matrix(ref);
    get_stl_string(str) = ss.str();
    ) 
}

extern "C"
int casadi_sx_matrix_print_cout(sx_matrix_ref ref){
  try{
    cout << get_sx_matrix(ref) << endl;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_print_cerr(sx_matrix_ref ref){
  try{
    cout << get_sx_matrix(ref) << endl;
    return 0;
  } catch(...){
    return 1;
  }
}

// extern "C"
// int casadi_sx_matrix_strbuf_len(){
//   return s.size();
// }
// 
// extern "C"
// const char* casadi_sx_matrix_strbuf(){
//   return s.c_str();
// }

extern "C"
int casadi_sx_matrix_binary(sx_matrix_ref r, int op, sx_matrix_ref x, sx_matrix_ref y){
  try{
    get_sx_matrix(r).binary(op,get_sx_matrix(x),get_sx_matrix(y));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_unary(sx_matrix_ref r, int op, sx_matrix_ref x){
  try{
    get_sx_matrix(r).unary(op,get_sx_matrix(x));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_prod(sx_matrix_ref r, sx_matrix_ref x, sx_matrix_ref y){
  try{
    get_sx_matrix(r) = prod(get_sx_matrix(x),get_sx_matrix(y));
    return 0;
  } catch(...){
    return 1;
  }    
}

extern "C"
int casadi_sx_matrix_vertcat(sx_matrix_ref r, sx_matrix_vec v){
  try{
    get_sx_matrix(r) = vertcat(get_sx_matrix_vec(v));    
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_horzcat(sx_matrix_ref r, sx_matrix_vec v){
  try{
    get_sx_matrix(r) = horzcat(get_sx_matrix_vec(v));    
    return 0;
  } catch(...){
    return 1;
  }
}


extern "C"
sx_matrix_vec casadi_sx_matrix_vec_new(void){
  try{
    return new vector<SXMatrix>();
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_sx_matrix_vec_delete(sx_matrix_vec v){
  try{
    vector<SXMatrix> *ptr = &get_sx_matrix_vec(v);
    delete ptr;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_vec_push_back(sx_matrix_vec v, sx_matrix_ref ref){
  try{
    get_sx_matrix_vec(v).push_back(get_sx_matrix(ref));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_vec_size(sx_matrix_vec v){
  try{
    return get_sx_matrix_vec(v).size();
  } catch(...){
    return -1;
  }
}

extern "C"
int casadi_sx_matrix_vec_print_cout(sx_matrix_vec v){
  try{
    cout << get_sx_matrix_vec(v) << endl;
    return 0;
  } catch(...){
    return 1;
  }
}


extern "C"
int casadi_sx_matrix_size(sx_matrix_ref ref, int *sz){
  try{
    if(sz==0) throw "null";
    *sz = get_sx_matrix(ref).size();
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_size1(sx_matrix_ref ref, int *sz){
  try{
    if(sz==0) throw "null";
    *sz = get_sx_matrix(ref).size1();
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_size2(sx_matrix_ref ref, int *sz){
  try{
    if(sz==0) throw "null";
    *sz = get_sx_matrix(ref).size2();
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_transpose(sx_matrix_ref res, sx_matrix_ref ref){
  try{
    get_sx_matrix(res) = trans(get_sx_matrix(ref));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_matrix_vec_print(sx_matrix_vec ref, string_ptr str){
  c_fcn_wrapper("casadi_sx_matrix_vec_print",
    stringstream ss;
    ss << get_sx_matrix_vec(ref);
    get_stl_string(str) = ss.str();
    )
}

extern "C"
int casadi_sx_matrix_sx(sx_matrix_ref ref, sx_ptr scalar){
  c_fcn_wrapper("casadi_sx_matrix_sx",
    get_sx_matrix(ref) = sx_ref(scalar);
  )
}

extern "C"
int casadi_sx_matrix_sx_vec(sx_matrix_ref ref, sx_vec_ptr v){
  c_fcn_wrapper("casadi_sx_matrix_sx_vec",
    get_sx_matrix(ref) = SXMatrix(sx_vec(v));
  )
}



} // namespace CasADi

