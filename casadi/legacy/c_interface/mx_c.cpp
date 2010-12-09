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

#include "mx_c.hpp"
#include "stl_string_c.hpp"
#include "../fx/mx_function.hpp"
#include "../stl_vector_tools.hpp"

#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>

using namespace CasADi;
using namespace std;

// Get a reference
MX& get_mx(mx_ref ref){
  if(ref==0) throw CasadiException("get_mx failed: null pointer");
  MX *r = (MX*)(ref);
  return *r;
}

// Get a reference to a vector
vector<MX>& get_mx_vec(mx_vec v){
  if(v==0) throw CasadiException("get_mx_vec failed: null pointer");
  vector<MX> *r = (vector<MX>*)(v);
  return *r;
}

extern "C"
mx_ref casadi_mx_new(){
  return new MX();
}

extern "C"
int casadi_mx_symbol(mx_ref ref, const char* name, int nrow, int ncol){
  try{
    get_mx(ref) = MX(name,nrow,ncol);
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_constant(mx_ref ref, const double* data, int nrow, int ncol, char order){
  try{
    get_mx(ref) = MX(data,nrow,ncol,order);
    return 0;
  } catch(...){
    return 1;
  }  
}

extern "C"
int casadi_mx_delete(mx_ref ref){
  MX* ptr = &get_mx(ref);
  delete ptr;
  return 0;
}

extern "C"
int casadi_mx_print_string(mx_ref ref, string_ptr str){
  try{
    // Print to a temporary variable
    stringstream ss;
    ss << get_mx(ref);
    get_stl_string(str) = ss.str();
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_print_cout(mx_ref ref){
  try{
    cout << get_mx(ref) << endl;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_print_cerr(mx_ref ref){
  try{
    cout << get_mx(ref) << endl;
    return 0;
  } catch(...){
    return 1;
  }
}

// extern "C"
// int casadi_mx_strbuf_len(){
//   return s.size();
// }
// 
// extern "C"
// const char* casadi_mx_strbuf(){
//   return s.c_str();
// }

extern "C"
int casadi_mx_binary(mx_ref r, int op, mx_ref x, mx_ref y){
  try{
    get_mx(r) = MX::binary(op,get_mx(x),get_mx(y));
    return 0;
  } catch(...){
    return 1;
  }    
}

extern "C"
int casadi_mx_unary(mx_ref r, int op, mx_ref x){
  try{
    get_mx(r) = MX::unary(op,get_mx(x));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_prod(mx_ref r, mx_ref x, mx_ref y){
  try{
    get_mx(r) = prod(get_mx(x),get_mx(y));
    return 0;
  } catch(...){
    return 1;
  }    
}

extern "C"
int casadi_mx_vertcat(mx_ref r, mx_vec v){
  try{
    get_mx(r) = vertcat(get_mx_vec(v));    
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_horzcat(mx_ref r, mx_vec v){
  try{
    get_mx(r) = horzcat(get_mx_vec(v));    
    return 0;
  } catch(...){
    return 1;
  }
}


extern "C"
mx_vec casadi_mx_vec_new(void){
  try{
    return new vector<MX>();
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_mx_vec_delete(mx_vec v){
  try{
    vector<MX> *ptr = &get_mx_vec(v);
    delete ptr;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_vec_push_back(mx_vec v, mx_ref ref){
  try{
    get_mx_vec(v).push_back(get_mx(ref));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_vec_size(mx_vec v){
  try{
    return get_mx_vec(v).size();
  } catch(...){
    return -1;
  }
}

extern "C"
int casadi_mx_vec_print_cout(mx_ref v){
  try{
    cout << get_mx_vec(v) << endl;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_size(mx_ref ref, int *sz){
  try{
    if(sz==0) throw "null";
    *sz = get_mx(ref).size();
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_size1(mx_ref ref, int *sz){
  try{
    if(sz==0) throw "null";
    *sz = get_mx(ref).size1();
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_size2(mx_ref ref, int *sz){
  try{
    if(sz==0) throw "null";
    *sz = get_mx(ref).size2();
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_mx_transpose(mx_ref res, mx_ref ref){
  try{
    get_mx(res) = trans(get_mx(ref));
    return 0;
  } catch(...){
    return 1;
  }
}



