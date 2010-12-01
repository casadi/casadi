#include "convection_diffusion_parameters.hpp"
#include "convection_diffusion_model.hpp"
#include <fstream>
#include <casadi/stl_vector_tools.hpp>

using namespace std;
using namespace OPTICON;


int main(){
  try{

  // create a model instance
  CDModel mod(NZ);

  // Vector for storing the temperature profile
  vector<double> t_prof;

  // Create a file for saving the results
  ofstream resfile;

  // Simulate the model and save the results to disk
  resfile.open ("results_convection_diffusion.txt");
  mod.simulate(resfile, t_prof);
  resfile.close();

  // save profile
  resfile.open ("t_prof.txt");
  write_matlab(resfile,t_prof);
  resfile.close();
  
  return 0;
} catch (const char * str){
  cerr << str << endl;
  return 1;
}

}
