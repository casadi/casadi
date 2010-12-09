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

#include "convection_diffusion_parameters.hpp"
#include "convection_diffusion_model.hpp"
#include <fstream>
#include <casadi/stl_vector_tools.hpp>

using namespace std;
using namespace OPTICON;


int main(){
  try{

  // read tabulated data from disk
  ifstream indata;

  // Read optimal trajectory from disk   
  vector<vector<double> > u_opt;
  indata.open ("optimal_u.txt");
  read_matlab(indata,u_opt);
  indata.close();

/*  vector<vector<double> > x_opt;
  indata.open ("optimal_x.txt");
  read_matlab(indata,x_opt);
  indata.close();*/
    
  // create a model instance
  CDModel mod(NZ);

  // create integrators
  Integrator_old integrator_pos = mod.make_integrator(1);
  Integrator_old integrator_neg = mod.make_integrator(-1);
  Integrator_old integrator_zero = mod.make_integrator(0);  

  // Create a file for saving the results
  ofstream resfile;
  resfile.open ("results_convection_diffusion_optimal.txt");

  // Vector for storing the temperature profile
  vector<double> t_prof(NZ,0.4);

  for(int i=0; i<u_opt.size(); ++i){
    double u = u_opt[i][0];

    OCP_old &ivp = u>0 ? mod.ivp_pos : u<0 ? mod.ivp_neg : mod.ivp_zero;
    Integrator_old &integrator = u>0 ? integrator_pos : u<0 ? integrator_neg : integrator_zero;
    
    // Pass the initial condition to the integrator
    const vector<int> &ic_ind = ivp.getVarInds(mod.T0);
    ivp.setArg(ic_ind,&t_prof[0]);
    
    // Pass the control to the integrator
    const vector<int> &u_ind = ivp.getVarInds(mod.u);
    ivp.setArg(u_ind,&u);

    // Integrate
    integrator->evaluate();

    // Get the results
    t_prof = integrator->output[0].data();
  
    // Save to disk
    for(vector<double>::const_iterator it=t_prof.begin(); it!=t_prof.end(); ++it)
      resfile << *it << " ";
    resfile << endl;
  }
  
   resfile.close();
  
  

  
  
  
// 
//       // Pass the control action
//       ivp[iu]->setArg(u_ind[iu], &u_disc[iu]);
//       
//       // Integrate
//       try{
//         integrator[iu]->evaluate();
//       } catch(...){
//          // The integration failed: set the value to be not a number
//           for(int l=0; l<n_dim; ++l){
// 	    assert(0);
// //            transition[iss][iu][l] = numeric_limits<double>::quiet_NaN();
// 	  }
//       }
// 
//       // New point in the state space
//       vector<double> x_new = x_cur;
//       
//       // doesn't work well, give some values:
//       x_new[0] = 0.5;
//       x_new[1] = 0.5;
// 
//       // Temperature profile
//       const vector<double> &T_prof_new = integrator[iu]->getOutput();

  
  
  
  
  return 0;
} catch (exception &e){
  cerr << e.what() << endl;  
  return 1;
} catch (const char * str){
  cerr << str << endl;
  return 1;
}

}
