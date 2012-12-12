#include <symbolic/casadi.hpp>
#include <iostream>

using namespace CasADi;
using namespace std;

int main(){
  int nEnd = 1000;
  for(int dIn=1; dIn<nEnd; ++dIn){
    cout << "trying "<< dIn << "/" << nEnd << "inputs: " << (100.0*dIn)/nEnd << "\%" << endl;
    SXMatrix xIn = ssym("x",dIn);
    FX f = SXFunction(xIn,pow(mul(trans(xIn),SXMatrix::ones(dIn)),2.));
    f.init();
    FX h = f.hessian();
  }
  
  return 0;

}