/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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
#include "base.hpp"
#include <stdint.h>
#include <string>
#include <cstring>
#include <typeinfo>



extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  if(nlhs!=0){
    mexErrMsgTxt("Expected nlhs = 0.");
  }

  if(nrhs!=1){
    mexErrMsgTxt("Expected nrhs = 1.");
  }
  
  const mxArray* in = prhs[0];
  if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS || mxIsComplex(in))
    mexErrMsgTxt("Input must be a real uint64 scalar.");
  Base* b = reinterpret_cast<Base*>(*((uint64_t *)mxGetData(in)));
  
  int a = b->foo();
  mexPrintf("a = %d \n",a);
}

