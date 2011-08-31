#include "matlab_interface.h"
#include <iostream>
#include <string>

using namespace std;

double test(double aa){
  cout << aa << endl;
  return aa*aa;
}

void* convert_to_swig(mxArray *array_ptr){
  cout << "convert_to_swig: " << array_ptr << endl;

  
  cout << "mxUNKNOWN_CLASS " << mxUNKNOWN_CLASS << endl;
  cout << "mxCELL_CLASS " << mxCELL_CLASS << endl;
  cout << "mxSTRUCT_CLASS " << mxSTRUCT_CLASS << endl;
  cout << "mxLOGICAL_CLASS " << mxLOGICAL_CLASS << endl;
  cout << "mxCHAR_CLASS " << mxCHAR_CLASS  << endl;
  cout << "mxDOUBLE_CLASS " << mxDOUBLE_CLASS  << endl;
  cout << "mxSINGLE_CLASS " << mxSINGLE_CLASS  << endl;
  cout << "mxINT8_CLASS " << mxINT8_CLASS << endl;
  cout << "mxUINT8_CLASS " << mxUINT8_CLASS << endl;
  cout << "mxSINGLE_CLASS " << mxSINGLE_CLASS  << endl;
  cout << "mxINT16_CLASS " << mxINT16_CLASS << endl;
  cout << "mxUINT16_CLASS " << mxUINT16_CLASS << endl;
  cout << "mxINT32_CLASS " << mxINT32_CLASS << endl;
  cout << "mxUINT32_CLASS " << mxUINT32_CLASS << endl;
  cout << "mxINT64_CLASS " << mxINT64_CLASS << endl;
  cout << "mxUINT64_CLASS " << mxUINT64_CLASS << endl;
  cout << "mxFUNCTION_CLASS" << mxFUNCTION_CLASS<< endl;
  cout << endl;

  cout << "mxGetClassID = " << mxGetClassID(array_ptr) << endl;
  cout << "mxGetClassName = " << mxGetClassName(array_ptr) << endl;

  string className = mxGetClassName(array_ptr);
  if(className.compare("swig_ref")==0){
    cout << "is a swig_ref pointer" << endl;
    int nlhs=0;
    mxArray *plhs[0] = {};
    int nrhs=1;
    mxArray *prhs[1] = {array_ptr};
    int flag = mexCallMATLAB(nlhs,plhs,nrhs,prhs,"convert_input_swig");
  }
  
  
  cout << mxGetM(array_ptr) << " x " << mxGetN(array_ptr) << endl;
  
  return const_cast<mxArray *>(array_ptr);
}

mxArray* convert_from_swig(void* proxy_ptr){
  cout << "convert_from_swig: " << proxy_ptr << endl;
  mxArray* ret = mxCreateDoubleMatrix(4, 2, mxREAL);
  return ret;
}

void convert_input_mx(void* proxy_ptr){
  cout << "convert_input_mx: " << proxy_ptr << endl;
}

void convert_input_swig(void* proxy_ptr){
  cout << "convert_input_swig: " << proxy_ptr << endl;
}



// void* convert_input2(const mxArray *array_ptr){
//   return const_cast<mxArray *>(array_ptr);
// }