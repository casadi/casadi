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

#include "stl_string_c.hpp"

using namespace std;

std::string& get_stl_string(string_ptr str){
  if(str==0) throw "getmxvec failed: null pointer";
  string *r = (string*)(str);
  return *r;  
}

extern "C"
string_ptr casadi_string_new(void){
  try{
    return new string();
  } catch(...){
    return 0;
  }
}

extern "C"
int casadi_string_delete(string_ptr str){
  try{
    string *ptr = &get_stl_string(str);
    delete ptr;
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_string_assign(string_ptr str, const char* s){
  try{
    get_stl_string(str) = string(s);
  } catch(...){
    return 1;
  }
}

extern "C"
const char* casadi_string_get(string_ptr str){
  try{
    return get_stl_string(str).c_str();
  } catch(...){
    return 0;
  }
}


