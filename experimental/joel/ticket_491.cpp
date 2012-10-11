#include <symbolic/casadi.hpp>

using namespace CasADi;
using namespace std;

int main(){
    
  SXMatrix u = ssym("u");
  SXMatrix x = ssym("x",2);
  SXMatrix f = x + u;
  f[0] += 1/x[1];
  f[1] += x[0];

  // Integrator
  vector<SXMatrix> ux(2);
  ux[0] = u;
  ux[1] = x;
  SXFunction F(ux,f);
  F.setOption("name","F");
  F.init();
  
  MX U("U");
  vector<double> X0(2);
  X0[0] = 0;
  X0[1] = 1;
  vector<MX> UX(2);
  UX[0] = U;
  UX[1] = X0;
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
//     gfcn.print();
    FX J = gfcn.jacobian();
    J.init();
 //   J.print();
    J.setInput(1);
    J.evaluate();
    cout << J.output() << endl;
    
  }
  
  return 0;

}