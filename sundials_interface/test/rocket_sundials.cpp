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

#include "sundials_interface/cvodes_integrator.hpp"
#include "sundials_interface/idas_integrator.hpp"
#include <casadi/stl_vector_tools.hpp>
#include <casadi/fx/simulator.hpp>
#include <casadi/fx/c_function.hpp>
#include "sundials_interface/cvodes_internal.hpp"

#include <fstream>
#include <iostream>

using namespace CasADi;
using namespace std;

// Use CVodes
const bool implicit_integrator = true;

// use plain c instead of SX
const bool plain_c = true;

// test adjoint sensitivities
const bool with_asens = false;

// The DAE residual in plain c
void dae_res_c(double tt, const double *yy, const double* yydot, const double* pp, double* rhs){
  // Get the arguments
  double s = yy[0], v = yy[1], m = yy[2];
  double u = pp[0];
  double sdot = yydot[0], vdot = yydot[1], mdot = yydot[2];

  // Calculate the DAE residual
  rhs[0] = sdot - v;
  rhs[1] = vdot - (u-0.02*v*v)/m;
  rhs[2] = mdot - (-0.01*u*u);
}

// Wrap the function to allow creating an CasADi function
void dae_res_c_wrapper(CFunction &f, int fsens_order, int asens_order, void* user_data){
  dae_res_c(f.input(DAE_T).data()[0],
            &f.input(DAE_Y).data()[0],
            &f.input(DAE_YDOT).data()[0],
            &f.input(DAE_P).data()[0],
            &f.output(DAE_RES).data()[0]);
}

int main(){
  // Time horizon
  double t0 = 0,  tf = 10;
  
  // Time 
  SX t("t");

  // Differential states
  SX s("s"), v("v"), m("m");
  vector<SX> y(3); 
  y[0] = s;
  y[1] = v;
  y[2] = m;

  // Control
  SX u("u");
  
  // Bounds on the control
  double u_lb = -0.5, u_ub = 1.3, u_init = 1;

  // Differential equation
  vector<SX> rhs(3);
  rhs[0] = v;              // sdot
  rhs[1] = (u-0.02*v*v)/m; // vdot
  rhs[2] = -0.01*u*u;      // mdot

  // Initial conditions
  vector<double> y0(3);
  y0[0] = 0;
  y0[1] = 0;
  y0[2] = 1;
  
  // Reference trajectory
  SX u_ref = 3-sin(t);
  
  // Square deviation from the state trajectory
  SX u_dev = u-u_ref;
  u_dev *= u_dev;

  // Full state (includes quadratures)
  vector<double> x0=y0;
  x0.push_back(0);

  // Integrator
  Integrator integrator;
  
  if(implicit_integrator){
    // Implicit integrator (IDAS)
    
    // State derivative
    SX sdot("sdot"), vdot("vdot"), mdot("mdot");
    vector<SX> ydot(3); 
    ydot[0] = sdot;
    ydot[1] = vdot;
    ydot[2] = mdot;
    
    // Input of the dae residual
    vector<vector<SX> > ffcn_in(DAE_NUM_IN);
    ffcn_in[DAE_T] = vector<SX>(1,t);
    ffcn_in[DAE_Y] = y;
    ffcn_in[DAE_YDOT] = ydot;
    ffcn_in[DAE_P] = vector<SX>(1,u);

    // DAE residual function
    FX ffcn;
    
    if(plain_c){
      // Use DAE residual defined in a c-function
      ffcn = CFunction(dae_res_c_wrapper);
      
      // Specify the number of inputs and outputs
      ffcn.setNumInputs(DAE_NUM_IN);
      ffcn.setNumOutputs(DAE_NUM_OUT);
      
      // Specify dimensions of inputs and outputs
      ffcn.input(DAE_T).setSize(1);
      ffcn.input(DAE_Y).setSize(3);
      ffcn.input(DAE_YDOT).setSize(3);
      ffcn.input(DAE_P).setSize(1);
      ffcn.output(DAE_RES).setSize(3);

    } else {
      // Use DAE residual evalable symbolically
      std::vector<SX> res = ydot;
      for(int i=0; i<res.size(); ++i) res[i] -= rhs[i];
      ffcn = SXFunction(ffcn_in,res);
      ffcn.setOption("ad_order",1);
    }
    ffcn.setOption("name","ODE right hand side");
    
    // Quadrature function
    SXFunction qfcn(ffcn_in,u_dev);
    qfcn.setOption("ad_order",1);

    // Create an integrator
    integrator = Sundials::IdasIntegrator(ffcn,qfcn);
    
    integrator.setOption("calc_ic",false);
    
  } else {
    // Explicit integrator (CVODES)
    
    // Input of the ode rhs
    vector<vector<SX> > ffcn_in(ODE_NUM_IN);
    ffcn_in[ODE_T] = vector<SX>(1,t);
    ffcn_in[ODE_Y] = y;
    ffcn_in[ODE_P] = vector<SX>(1,u);
  
    // ODE right hand side
    SXFunction ffcn(ffcn_in,rhs);
    ffcn.setOption("name","ODE right hand side");
    ffcn.setOption("ad_order",1);
    
    // Quadrature function
    SXFunction qfcn(ffcn_in,u_dev);
    qfcn.setOption("ad_order",1);
  
    // Create an integrator
    integrator = Sundials::CVodesIntegrator(ffcn,qfcn);
  }

  integrator.setOption("ad_order",1);
  integrator.setOption("fsens_err_con",true);
  integrator.setOption("quad_err_con",true);
  integrator.setOption("abstol",1e-12);
  integrator.setOption("reltol",1e-12);
  integrator.setOption("fsens_abstol",1e-6);
  integrator.setOption("fsens_reltol",1e-6);
  integrator.setOption("asens_abstol",1e-6);
  integrator.setOption("asens_reltol",1e-6);
  if(plain_c){
    //integrator.setOption("exact_jacobian",false);
    integrator.setOption("finite_difference_fsens",true);
  } else {
    //integrator.setOption("exact_jacobian",true);
    integrator.setOption("finite_difference_fsens",false);    
  }

  integrator.setOption("max_num_steps",100000);
  //  integrator.setOption("linear_solver","iterative");
  
  //  cout << "Integrator:" << endl;
//  integrator.printOptions();

  if(1){
  integrator.init();
  
//   integrator.input(INTEGRATOR_T0).set(t0);
//   integrator.input(INTEGRATOR_TF).set(tf);
//   integrator.input(INTEGRATOR_X0).set(x0);
//   integrator.input(INTEGRATOR_P).set(u_init);

  integrator.setInput(t0,INTEGRATOR_T0);
  integrator.setInput(tf,INTEGRATOR_TF);
  integrator.setInput(x0,INTEGRATOR_X0);
  integrator.setInput(u_init,INTEGRATOR_P);
  
  double yp0[] = {0,1,-0.01,0};
  integrator.setInput(yp0,INTEGRATOR_XP0);
  
  integrator.evaluate();
  cout << "before " << integrator.getOutputData() << endl;
  vector<double> res0 = integrator.getOutputData();

//  x0[1] += 0.01;
  u_init += 0.01;

  integrator.setInput(x0,INTEGRATOR_X0);
  integrator.setInput(u_init,INTEGRATOR_P);
  integrator.evaluate();
  
  // Print statistics
  integrator.printStats();

  cout << "after " << integrator.output().data() << endl;
  
  vector<double> fd = integrator.output().data();
  for(int i=0; i<fd.size(); ++i){
    fd[i] -= res0[i];
    fd[i] /= 0.01;
  }
  
  cout << "fd    " << fd << endl;

  vector<double> x0_seed(x0.size(),0);
  double u_seed = 0;
  
  u_seed = 1;
//  x0_seed[1] = 1;

  
  integrator.setInput(u_init,INTEGRATOR_P);

  // forward seeds
  integrator.setFwdSeed(0.0,INTEGRATOR_T0);
  integrator.setFwdSeed(0.0,INTEGRATOR_TF);
  integrator.setFwdSeed(x0_seed,INTEGRATOR_X0);
  integrator.setFwdSeed(u_seed,INTEGRATOR_P);

  if(with_asens){
    // backward seeds
    vector<double> &bseed = integrator.output(INTEGRATOR_XF).dataA();
    fill(bseed.begin(),bseed.end(),0);
    bseed[0] = 1;

    // evaluate with forward and adjoint sensitivities
    integrator.evaluate(1,1);
  } else {
    // evaluate with only forward sensitivities
    integrator.evaluate(1,0);
  }
    
  vector<double> &fsens = integrator.output().dataF();
  cout << "fsens " << fsens << endl;

  if(with_asens){
    cout << integrator.input(INTEGRATOR_T0).dataA() << endl;
    cout << integrator.input(INTEGRATOR_TF).dataA() << endl;
    cout << integrator.input(INTEGRATOR_X0).dataA() << endl;
    cout << integrator.input(INTEGRATOR_P).dataA() << endl;
  }
  
  return 0;
  }
  
//  vector<double> vv(1);
//  setv(0,vv);
  
  // Create a grid
  int ntimes = 20;
  // The desired output times
  vector<double> t_out(ntimes);
  linspace(t_out,0,10);
  
  // Create a simulator
  Simulator simulator(integrator,t_out);
  simulator.setOption("name","rocket simulator");
  simulator.setOption("ad_order",1);
  
  cout << "Simulator:" << endl;
  simulator.printOptions();
  
  // Initialize the simulator to allow evaluation
  simulator.init();
  
  // Pass initial condition
  simulator.setInput(x0,SIMULATOR_X0);

  // Pass parameter values
  simulator.setInput(u_init,SIMULATOR_P);
  
  // Integrate
  simulator.evaluate();
  
  // Print to screen
  vector<double> xout(x0.size()*ntimes);
  simulator.getOutput(xout);
  cout << "xout = " << xout << endl;

  // Finite difference approximation
  simulator.setInput(u_init+0.01,SIMULATOR_P);
  simulator.evaluate();
  vector<double> fd(xout.size());
  simulator.getOutput(fd);
  for(int i=0; i<xout.size(); ++i){
    fd[i] -= xout[i];
    fd[i] /= 0.01;
  }
  cout << endl << "finite difference approximation =   " << fd << endl;
  
  // Calculate internally with forward sensitivities
  simulator.setFwdSeed(1.0,SIMULATOR_P);
  
  // Integrate with forward sensitivities
  simulator.evaluate(1,0);
  
  // Print to screen
  vector<double> xsens(x0.size()*ntimes);
  simulator.getFwdSens(xsens);
  cout << endl << "forward sensitivity result = " << xsens << endl;

  // Pass backward seeds
  vector<double> bsens(xsens.size(),0.0);
  bsens[x0.size()*(ntimes-1)+1] = 1;
  simulator.setAdjSeed(bsens);

  // Integrate backwards
  simulator.evaluate(0,1);
  
  // Print results
  cout << "adjoint sensitivity result (p): " << simulator.input(SIMULATOR_P).data() << endl;
  cout << "adjoint sensitivity result (x0): "<< simulator.input(SIMULATOR_X0).data() << endl;
  
  
  #if 0
  // Parametrize controls into 20 uniform intervals
  Matrix u_disc = parametrizeControls(ocp,u,20);
//  Matrix udisc_guess = 2*ones(20,1);
//  ocp.guessSolution(u_disc, udisc_guess);


  // numeric
  vector<double> t_out_num(nt_out);
  for(int i=0; i<nt_out; ++i)
    t_out_num[i] = i*10/double(nt_out-1);
 
  // LAGRANGE TERM
//  int l_ind = ocp.addOutput(t_out, Matrix(), u*u, 1, u_disc); // quadrature and sensitivities, first order
//  int l_ind = ocp.addOutput(10, Matrix(), u*u); // just quadrature 
//  int l_ind = ocp.addOutput(10, 100*(1-m),Matrix()); // just quadrature 
  int l_ind = ocp.addOutput(10, m, Matrix(), 1, u_disc); // quadrature and sensitivities, first order

  // forward sensitivities of v with respect to u at time tf
//   ocp.addOutput(10, v, Matrix(), 1, u_disc);

  
  // Set output function
  int oind_s = ocp.addOutput(t_out, s);
  int oind_v = ocp.addOutput(t_out, v);
  int oind_m = ocp.addOutput(t_out, m);
  int oind_u = ocp.addOutput(t_out, u);

  // Eliminate dependent variables from the functions
  eliminateDependent(ocp);

  // Print to screen
  std::cout << ocp;

  // Allocate an integrator
  CvodesIntegrator integrator(ocp);

  // Print the possible options
  integrator->printOptions();

  // Set the linear solver
   integrator->setOption("linear_solver","dense");
//  integrator->setOption("linear_solver","band");
//  integrator->setOption("linear_solver","sparse");

  // Use exact jacobian
  integrator->setOption("exact_jacobian",false);

  // Upper and lower band-widths (only relevant when using a band linear solver)
//   integrator->setOption("mupper",1);
//   integrator->setOption("mlower",1);

  // set tolerances 
  integrator->setOption("reltol",1e-6);
  integrator->setOption("abstol",1e-8);

  // Initialize the integrator
  integrator->init();
  std::cout << "initialized" << std::endl;

  // Integrate once for visualization
  integrator->evaluate();
  std::cout << "solved" << std::endl;

  // Create a file for saving the results
  std::ofstream resfile;
  resfile.open ("results_rocket.txt");



  // Save results to file
  resfile << "t_out " << t_out_num << std::endl;
  resfile << "s " << integrator->output[oind_s].data() << endl;
  resfile << "v " << integrator->output[oind_v].data() << endl;
  resfile << "m " << integrator->output[oind_m].data() << endl;
  resfile << "u " << integrator->output[oind_u].data() << endl;
  resfile << "lfun " << integrator->output[l_ind].data() << endl;
  cout    << "lfun " << integrator->output[l_ind].data() << endl;

  // Close the results file
  resfile.close();

  return 0;
} catch (const char * str){
  std::cerr << str << std::endl;
  return 1;
}
#endif

}
