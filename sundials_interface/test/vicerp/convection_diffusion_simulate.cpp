#include <fstream>
#include <ctime>
#include "convection_diffusion_model.hpp"




int main(){
try {

   clock_t time0 = clock();
   cout << "program started " << endl;

  // Calculate the derivative approximations
  calc_Tdisc();

  // Create a file for saving the results
  ofstream resfile;
  resfile.open ("results_convection_diffusion.txt");

  // Vector for storing the temperature profile
  vector<double> t_prof;

  // Simulate the model
  simulate(resfile, t_prof);

  return 0;
} catch (const char * str) {
  cerr << str << endl;
  return 1;
} catch (exception &e) {
  cerr << "fatal error: " << e.what() << endl;
  return 1;
}

}



