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

#include "mex.h"
#include "ipopt_solver_c.h"
#include "ipopt_solver.hpp"
#include "../casadi/c_interface/fx_c.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

using namespace CasADi;
using namespace std;

extern "C"
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
 if (nlhs != 0 || nrhs != 5)
 mexErrMsgTxt("Wrong number of arguments.");
 
 /* Get the pointers */
 int *fcn = (int*)mxGetData(prhs[0]);
 int *ffcn = (int*)mxGetData(prhs[1]);
 int *gfcn = (int*)mxGetData(prhs[2]);
 int *hfcn = (int*)mxGetData(prhs[3]);
 int *jfcn = (int*)mxGetData(prhs[4]);
  
  // this is never reached
  int flag = casadi_ipopt_solver((void*)fcn[0], (void*)ffcn[0], (void*)gfcn[0], (void*)hfcn[0], (void*)jfcn[0]);
  if(flag!=0) 
    mexErrMsgTxt("Error in casadi_ipopt_solver");

}




