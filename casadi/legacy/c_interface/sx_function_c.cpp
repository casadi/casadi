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

#include "sx_function_c.h"
#include "../fx/sx_function.hpp"
#include "sx_matrix_c.hpp"
#include "fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

using namespace CasADi;
using namespace std;


extern "C"
int casadi_sx_function(fx_ref fcn, sx_matrix_vec iv, sx_matrix_vec ov){
 try{
    get_fx(fcn) = SXFunction(get_sx_matrix_vec(iv),get_sx_matrix_vec(ov));
    return 0;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_sx_function_evaluate_symbolically(fx_ref fcn, sx_matrix_vec iv, sx_matrix_ref res){
 try{
    SXFunction fun(get_fx(fcn));
    get_sx_matrix(res) = fun(get_sx_matrix_vec(iv));
    return 0;
  } catch(...){
    return 1;
  }
}
