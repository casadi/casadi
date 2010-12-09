/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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



