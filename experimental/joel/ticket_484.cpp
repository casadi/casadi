#include <symbolic/casadi.hpp>

using namespace CasADi;
using namespace std;

int main(){
  MX x("x");
  MXFunction f(x,x);
  f.setOption("ad_mode","reverse");
  f.setOption("name","f");
  f.init();
  FX J = f.jacobian();
  f.init();
  J.init();
  J.evaluate();
  
  return 0;
}