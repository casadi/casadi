#include <symbolic/casadi.hpp>

using namespace CasADi;
using namespace std;

int main(){
    
  // Create simple function
  vector<SXMatrix> ux(2);
  SXMatrix u = ux[0] = ssym("u");
  SXMatrix x = ux[1] = ssym("x");
  SXFunction F(ux,u + 1/x);
  F.setOption("name","F");
  F.init();
  
  // Call this function in a graph
  MX U("U");
  vector<MX> UX(2);
  UX[0] = U;
  UX[1] = U;
  MX X = F.call(UX).at(0);
  UX[0] = U;
  UX[1] = X;
  MX G = F.call(UX).at(0);

  for(int kk=0; kk<3; ++kk){
    switch(kk){
      case 0: cout << "SX SCT" << endl; break;
      case 1: cout << "MX OO" << endl; break;
      case 2: cout << "MX SCT" << endl; break;
    }
    
    FX gfcn;
    if(kk==0){
      MXFunction tmp(U,G);
      tmp.init();
      gfcn = tmp.expand();
    } else {
      gfcn = MXFunction(U,G);
      gfcn.setOption("numeric_jacobian",kk==1);
    }
    gfcn.setOption("ad_mode","reverse");
    
    gfcn.init();
    gfcn.print();
    FX J = gfcn.jacobian();
    J.init();
 //   J.print();
    J.setInput(1);
    J.evaluate();
    cout << J.output() << endl;
    
  }
  
  return 0;

}