#include <casadi/casadi.hpp>

using namespace CasADi;
using namespace std;

int main(){
  
  // Construct a simple function
  SXMatrix x1 = ssym("x1");
  SXMatrix x2 = ssym("x2");
  SXMatrix r1 = sin(x2);
  SXMatrix r2 = x1+5;
  
  // Input arguments
  vector<SXMatrix> F_in(2);
  F_in[0] = x1;
  F_in[1] = x2;
  
  // Output arguments
  vector<SXMatrix> F_out(2);
  F_out[0] = r1;
  F_out[1] = r2;
  
  // Create function
  SXFunction F(F_in,F_out);
  F.setOption("just_in_time",true);
  F.init();

  return 0;
} 
