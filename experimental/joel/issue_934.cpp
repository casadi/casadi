#include "symbolic/casadi.hpp"
#include "integration/irk_integrator.hpp"
#include "interfaces/sundials/kinsol_solver.hpp"

using namespace CasADi;
using namespace std;

// Print a (typically 64-bit) unsigned integer
void printBinary(bvec_t v){
  for(int k=0; k<bvec_size; ++k){
    if(k%4==0) cout << " ";
    if(v & (bvec_t(1)<<k)){
      cout << 1;
    } else {
      cout << 0;
    }
  }
  cout << endl;
}

int main(){
  // DAE
  SXMatrix x = ssym("x",3);
  SXMatrix ode = SXMatrix::zeros(3);
  ode[0] = x[0];
  ode[1] = x[0]+x[1];
  ode[2] = x[0]+x[1]+x[2];
  SXFunction dae(daeIn("x",x),daeOut("ode",ode));

  // Integrator
  IRKIntegrator f(dae);
  f.setOption("implicit_solver",KinsolSolver::creator);
  f.init();
    
  // Get arrays for the inputs and outputs, reinterpreting the vector of double as an array of unsigned integers
  bvec_t* f_in = get_bvec_t(f.input(INTEGRATOR_X0).data());
  bvec_t* f_out = get_bvec_t(f.output(INTEGRATOR_XF).data());
    
  // Propagate from input to output (forward mode)
  int fwd = true;
    
  // Make sure that the class is able to support the dependency propagation
  casadi_assert(f.spCanEvaluate(fwd));
    
  // Pass seeds
  fill(f_in,f_in,0);
  f_in[0] = bvec_t(1) << 0; // seed in direction 0
  f_in[1] = bvec_t(1) << 2; // seed in direction 2
  f_in[2] = (bvec_t(1) << 4) | (bvec_t(1) << 63); // seed in direction 4 and 63

  // Reset sensitivities
  f_out[0] = 0;

  cout << "forward seeds" << endl;
  printBinary(f_in[0]);
  printBinary(f_in[1]);
  printBinary(f_in[2]);
    
  // Propagate dependencies
  f.spInit(fwd);
  f.spEvaluate(fwd);

  // Print the result
  cout << "forward sensitivities" << endl;
  printBinary(f_out[0]);
  printBinary(f_out[1]);
  printBinary(f_out[2]);
    
  return 0;
  // Propagate from output to input (adjoint/reverse/backward mode)
  cout << "backward mode" << endl;
  fwd = false;

  // Make sure that the class is able to support the dependency propagation
  casadi_assert(f.spCanEvaluate(fwd));

  // Pass seeds
  f_out[0] = (bvec_t(1) << 5) | (bvec_t(1) << 6); // seed in direction 5 and 6
    
  // Reset sensitivities
  f_in[0] = 0;
  f_in[1] = 0;
  f_in[2] = 0;
    
  // Propagate dependencies
  f.spInit(fwd);
  f.spEvaluate(fwd);

  // Print the result
  printBinary(f_in[0]);
  printBinary(f_in[1]);
  printBinary(f_in[2]);
  
  return 0;
}

