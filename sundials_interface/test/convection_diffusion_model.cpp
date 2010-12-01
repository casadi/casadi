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

#include "convection_diffusion_model.hpp"
#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>

CDModel::CDModel(int NZ){
   
  // parameters
  T_l = (T_receiver+273)/(T_return+273); // inlet temperature during loading
  T_u = 1; // inlet temperature during unloading
  D = diff; // diffusion coefficient
  
  // discretization
  SX dz = h/NZ;
  SX dz2 = dz*dz;

  // parameters
  SX T_l = (T_receiver+273)/(T_return+273); // inlet temperature during loading
  SX T_u = 1; // inlet temperature during unloading
    
    // differential states
  T = Matrix("T",NZ);
  Tdot = Matrix("Tdot",NZ);
      
  dT_dz_pos = (T[0] - T_l)/dz;
  for(int i=1; i<NZ; ++i){
    Matrix Tder_i = (T[i] - T[i-1])/dz;
    dT_dz_pos << Tder_i;
    dT_dz_neg << Tder_i;
  }
  dT_dz_neg << (T_u - T[NZ-1])/dz;

  d2T_dz2 << (-T[0] + T[1])/dz2;
  for(int i=1; i<NZ-1; ++i)
    d2T_dz2 << (T[i-1]- 2*T[i] + T[i+1])/dz2;
 
  d2T_dz2 << (T[NZ-2]- T[NZ-1])/dz2;
  

  // Control (flow to storage)
  u = SX("u");

  // Maximal u
  u_max = SX(0.33); // fully loading/unloading takes 1.5 hours

  // Initial temperature profile
  T0 = Matrix("T0",NZ);

}

// Just simulate a scenario
void CDModel::simulate_new(ostream& resfile, vector<double>& t_prof){
    // time points
  double t_sunrise = 0;
  double t_startup = 1;
  double t_loading_start = 6-sqrt(7.2);
  double t_loading_end = 6+sqrt(7.2);
  double t_shutdown = 11;
  double t_sunset = 12;
  double t_next_sunrise = 24;
  
  // Time variable
  SX t("t");
  
  // Solar irradiation
  SX phiR = (t<t_sunset) * t*(12-t)/36; // quadratic function

  // Mass flow to boiler
  SX phiB = 0.8 * (( t >= t_startup) && ( t < t_shutdown));

  // Differential states
  SX eta = 1 - ((0.80 - phiB)*(0.80 - phiB))/((0.80 - 0.25)*(0.80 - 0.25)); // quadratic function with max at 80 %
  
  // Mass flow to storage (shadow member)
  SX u = phiR - phiB;

  // Differential equation
  Matrix Tdot = - u * if_else(u>=0, dT_dz_pos, dT_dz_neg) + D*d2T_dz2;

  // Initial value
  vector<double> T0(NZ,1);

  // Create an ode rhs
  vector<Matrix> f_input(3);
  f_input[0] = t;
  f_input[1] = T;
  f_input[2] = Matrix(); // no parameter as for now
  SXFunction f(f_input,Tdot);
  
  // Create an integrator
  CVodesInterface integrator(f);
  integrator.setOption("linear_solver","banded");
//  integrator.setOption("iterative_solver","gmres");
  integrator.setOption("upper_bandwidth",1);
  
  integrator.setOption("lower_bandwidth",1);
  integrator.setOption("max_num_steps",1000);
  integrator.init();
  
  // Pass initial value
  integrator.setInput(&t_sunrise,0);
  integrator.setInput(&t_startup,1);
  integrator.setInput(T0,2);

  // Integrate
  cout << "stating simulation" << endl;
  integrator.evaluate();
  cout << "simulation ended" << endl;
  
  // Get output
  vector<double> T1(NZ);
  integrator.getOutput(T1);

  cout << "T1 = " << T1 << endl;
  
  // Print statistics
  integrator.printStats();
  
  
  
}
  
void CDModel::simulate(ostream& resfile, vector<double>& t_prof){
  simulate_new(resfile,t_prof);
    
  return;
  
    // time points
  SX t_sunrise = 0;
  SX t_startup = 1;
  SX t_loading_start = 6-sqrt(7.2);
  SX t_loading_end = 6+sqrt(7.2);
  SX t_shutdown = 11;
  SX t_sunset = 12;
  SX t_next_sunrise = 24;

  // Create a single-stage optimal control problem
  OCP_old ocp;

  // Time variable
  Matrix t = ocp.getTime();

  // Specify the time horizon
  ocp.setInitialTime(t_sunrise);
  ocp.setFinalTime(t_next_sunrise);

  // differential states
  ocp.addState(T,Tdot);

  // Solar irradiation
  Matrix phiR = (t<t_sunset) * t*(12-t)/36; // quadratic function

  // Mass flow to boiler
  Matrix phiB = 0.8 * (( t >= t_startup) && ( t < t_shutdown));

  // Differential states
  Matrix eta = 1 - ((0.80 - phiB)*(0.80 - phiB))/((0.80 - 0.25)*(0.80 - 0.25)); // quadratic function with max at 80 %
  
  // Mass flow to storage (shadow member)
  Matrix u = phiR - phiB;

  // Differential equation
  ocp.addEquation(Tdot, - u * if_else(u>=0, dT_dz_pos, dT_dz_neg) + D*d2T_dz2);

  // Initial temperature in the storage
  Matrix T0 = T_u * ones(NZ,1);
  ocp.addInitialCondition(T, T0);

  // The desired output times
  Matrix t_out = linspace(0,10,nt_out);

  // numerically
  vector<double> t_out_num(nt_out);
  for(int i=0; i<nt_out; ++i)
    t_out_num[i] = 10*i/double(nt_out-1);

  // Project to these points
  Matrix T_proj(nz_out,NZ);
  for(int i=0; i<nz_out; ++i)
    T_proj(i,(NZ*i)/nz_out) = 1; // integer division!

  // Set output function
  int oind_T = ocp.addOutput(t_out, T_proj*T);
  int oind_u = ocp.addOutput(t_out, u);

  // Eliminate dependent variables from the functions
  eliminateDependent(ocp);

  // Allocate an integrator
  CvodesIntegrator integrator(ocp);

  // Set the linear solver
  //   integrator.setOption("linear_solver","dense");
  integrator->setOption("linear_solver","band");
  //  integrator->setOption("linear_solver","sparse");

  // Use exact jacobian
  integrator->setOption("exact_jacobian",false);

  // Upper and lower band-widths (only relevant when using a band linear solver)
   integrator->setOption("mupper",1);
   integrator->setOption("mlower",1);

  // set tolerances 
  integrator->setOption("reltol",1e-6);
  integrator->setOption("abstol",1e-8);

  // Initialize the integrator
  integrator->init();
  cout << "initialized" << endl;

  // Integrate once for visualization
  integrator->evaluate();
  cout << "solved" << endl;

  // Save results to file
  resfile << "t_out " << t_out_num << endl;
  resfile << "T " << integrator->output[oind_T].data() << endl;
  resfile << "u " << integrator->output[oind_u].data() << endl;
  resfile << "nz_out " << nz_out << endl;
  resfile << "nt_out " << nt_out << endl;

  vector<double> z_out(nz_out);
  for(int i=0; i<nz_out; ++i)
    z_out[i] = i/double(nz_out-1);

  resfile << "z_out " << z_out << endl;

  // Save to output
  t_prof = integrator->output[oind_T].data();
}

CVodesInterface CDModel::make_integrator_new(int dir){
  // Time variable
  SX t("t");
  
  // Control
  Matrix f;
  if(dir==1){
      f = -u*u_max*dT_dz_pos + D*d2T_dz2;
  } else if(dir==-1) {
      f = -u*u_max*dT_dz_neg + D*d2T_dz2;      
  } else {
      f = D*d2T_dz2;      
  }
  
  // Create the ode rhs function
  vector<Matrix> input(3);
  input[0] = t;
  input[1] = T;
  input[2] = u; 
  SXFunction ffcn(input,f);
  
  // Create an integrator
  CVodesInterface integrator(ffcn);
  integrator.setOption("linear_solver","banded");
//  integrator.setOption("iterative_solver","gmres");
  integrator.setOption("upper_bandwidth",1);
  integrator.setOption("lower_bandwidth",1);
  integrator.setOption("max_num_steps",1000);
  integrator.init();

  // Set the time horizon
  double t0 = 0, tf = DT;
  integrator.setInput(&t0,0);
  integrator.setInput(&tf,1);
  
  
  return integrator;
}
  
CvodesIntegrator CDModel::make_integrator(int dir){
    OCP_old &ivp = dir>0 ? ivp_pos : dir<0 ? ivp_neg : ivp_zero;
  
    // Time variable
    Matrix t = ivp.getTime();
    
    // Specify the time horizon
    ivp.setInitialTime(0);
    ivp.setFinalTime(DT);
    
    // differential states
    ivp.addState(T,Tdot);

    // control
    ivp.addParameter(u);

    // Initial temperature in the storage
    ivp.addParameter(T0);
    ivp.addInitialCondition(T, T0);

    // Set output function (final T)
    ivp.addOutput(DT, T);

    if(dir==1){
      ivp.addEquation(Tdot, -u*u_max*dT_dz_pos + D*d2T_dz2);
    } else if(dir==-1) {
      ivp.addEquation(Tdot, -u*u_max*dT_dz_neg + D*d2T_dz2);      
    } else {
      ivp.addEquation(Tdot, D*d2T_dz2);      
    }
    
    CvodesIntegrator integrator = CvodesIntegrator(ivp);
    
    // Set the linear solver
    integrator->setOption("linear_solver","band");

    // Upper and lower band-widths (only relevant when using a band linear solver)
    integrator->setOption("mupper",1);
    integrator->setOption("mlower",1);

    // Use exact jacobian
    integrator->setOption("exact_jacobian",false);

    // set tolerances 
    integrator->setOption("reltol",1e-6);
    integrator->setOption("abstol",1e-8);

    // Initialize the integrator
    integrator->init();
    cout << "initialized" << endl;

    // Set options
    integrator->setOption("use_events",  false);
    integrator->setOption("use_output",  true);
    integrator->setOption("use_fsens",   false);
    integrator->setOption("use_asens",   false);
    integrator->setOption("use_quad",    false);
    
    return integrator;
}



