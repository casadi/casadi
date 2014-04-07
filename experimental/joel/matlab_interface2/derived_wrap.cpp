#include "mex.h" 
#include "derived.hpp"
#include <stdint.h>
#include <string>
#include <cstring>
#include <typeinfo>


extern "C" 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  if(nlhs!=1){
    mexErrMsgTxt("Expected nlhs = 1.");
  }

  if(nrhs!=0){
    mexErrMsgTxt("Expected nrhs = 0.");
  }

  // Create derived function
  Derived* d = new Derived();
  
  mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((uint64_t *)mxGetData(out)) = reinterpret_cast<uint64_t>(d);
  plhs[0] = out;
}


