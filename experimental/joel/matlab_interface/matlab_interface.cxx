#include "matlab_interface.h"
#include <iostream>

using namespace std;

double test(double aa){
  cout << aa << endl;
  return aa*aa;
}

int test_mex(const mxArray *array_ptr){
  cout << mxGetM(array_ptr) << " x " << mxGetN(array_ptr) << endl;
  
  
  return 5;
}
