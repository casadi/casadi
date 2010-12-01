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




