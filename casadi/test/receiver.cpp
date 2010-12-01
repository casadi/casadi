#include "../matrix.hpp"
#include "../expression_library.hpp"
#include "../stl_vector_tools.hpp"
#include "../sx_function.hpp"
#include <cassert>

using namespace std;


int main(){

try{
  // Dimensions
  //  int nl = 4*4; // number of absorber modules
  int nk = 200; // number of absorber cups per absorber module

  // Constants
  Matrix h_in(nk,1); // incoming specific enthalpy
  double A = 1; // absorber cup area

  // Free variables
  Matrix P("P",nk);  // Density
  Matrix dp("dp");   // Pressure drop
  Matrix mdot("mdot",nk); // mass flow through each absorber cup
  
  // Objective function
  Matrix f = -ones(1,nk)*mdot;

  // Efficiency
  Matrix eta(nk,1); 
  eta = mdot + P; //= ... 5-th degree polynomial of mdot and P

  // specific enthalpy for each absorber cup
  Matrix h_abs(nk,1); 
  for(int k=0; k<nk; ++k)
    h_abs[k] = h_in[k] + eta[k]*P[k]*A/mdot[k];
 
  // Pressure drop equation for each cup
  Matrix def_dp(nk,1);
  for(int k=0; k<nk; ++k){
    Matrix k1 = sin(h_abs[k]); // = ... (function of h_abs[k])
    Matrix k2 = atan(h_abs[k]); // = ... (function of h_abs[k])
    def_dp[k] = k1*mdot[k] + k2*mdot[k]*mdot[k];
  }
    
  // Gradient of f with respect to mdot
  Matrix grad_f = 

  
  
  
  

  


  return 0;

} catch (const char * str){
  std::cerr << str << std::endl;
  return 1;
}


}
