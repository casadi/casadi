#include <symbolic/casadi.hpp>

using namespace casadi;
using namespace std;

int main(){
  
  // Construct a simple function
  SX x1 = ssym("x1");
  SX x2 = ssym("x2");
  SX r1 = sin(x2);
  SX r2 = x1+5;
  
  // Input arguments
  vector<SX> F_in(2);
  F_in[0] = x1;
  F_in[1] = x2;
  
  // Output arguments
  vector<SX> F_out(2);
  F_out[0] = r1;
  F_out[1] = r2;
  
  // Create function
  SXFunction F(F_in,F_out);
  F.setOption("just_in_time",true);
  F.init();

  // Generate C code
  F.generateCode("test.c");
  
  // Pass inputs
  F.setInput(10,0);
  F.setInput(20,1);
  
  // Evaluate
  F.evaluate();
  
  // Print results
  cout << F.output(0) << endl;
  cout << F.output(1) << endl;
  
  // Print the LLVM IR
  F.print();
  
  return 0;
} 
