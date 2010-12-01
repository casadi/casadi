/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

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
