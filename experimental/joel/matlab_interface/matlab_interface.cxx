#include "matlab_interface.h"
#include <iostream>

using namespace std;

double test(double aa){
  cout << aa << endl;
  return aa*aa;
}

void* convert_to_swig(const mxArray *array_ptr){
  cout << mxGetM(array_ptr) << " x " << mxGetN(array_ptr) << endl;
  
  return const_cast<mxArray *>(array_ptr);
}

mxArray* convert_from_swig(void* proxy_ptr){
  mxArray* ret = mxCreateDoubleMatrix(4, 2, mxREAL);
  return ret;
}




// void* convert_input2(const mxArray *array_ptr){
//   return const_cast<mxArray *>(array_ptr);
// }