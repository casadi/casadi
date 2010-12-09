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

#include <acado_interface/acado_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>

#include <fstream>
#include <iostream>

using namespace CasADi;
using namespace std;

int main(){
  // Time horizon
  double t0 = 0,  tf = 10;
  
  // Time 
  SX t("t");

  // Differential states
  SX s("s"), v("v"), m("m");
  
  // Lagrangian term
  SX l("l");

  // State vector (includes quadratures)
  vector<SX> x;
  x.push_back(s);
  x.push_back(v);
  x.push_back(m);
  x.push_back(l);

  // Control
  SX u("u");
  
  // Bounds on the control
  double u_lb = -0.5, u_ub = 1.3, u_init = 1;

  // Reference trajectory
  SX u_ref = 3-sin(t);
  
  // Square deviation from the state trajectory
  SX u_dev = u-u_ref;
  u_dev *= u_dev;

  // Differential equation
  SX sdot= v;
  SX vdot= (u-0.02*v*v)/m;
  SX mdot= -0.01*u*u;
  SX ldot= u_dev;
  vector<SX> xdot;
  xdot.push_back(sdot);
  xdot.push_back(vdot);
  xdot.push_back(mdot);
  xdot.push_back(ldot);

  // Initial conditions
  vector<double> x0;
  x0.push_back(0);
  x0.push_back(0);
  x0.push_back(1);
  x0.push_back(0);
  
  // Input of the ode rhs
  vector<vector<SX> > ffcn_in(3);
  ffcn_in[0] = vector<SX>(1,t);
  ffcn_in[1] = x;
  ffcn_in[2] = vector<SX>(1,u);
  
  // ODE right hand side
  SXFunction ffcn(ffcn_in,xdot);
  ffcn.setOption("name","ODE right hand side");
  ffcn.setOption("ad_order",1);

  // Create an integrator
  ACADOIntegrator integrator(ffcn);
  integrator.setOption("ad_order",2);
  integrator.setOption("print_level","low");
  integrator.setOption("abstol",1e-10);
  integrator.setOption("reltol",1e-10);
  integrator.setOption("exact_jacobian",true);

  // Print options
  integrator.printOptions();

  // Initialize the integrator
  integrator.init();
  
  // Pass arguments
  integrator.input(INTEGRATOR_T0).set(t0);
  integrator.input(INTEGRATOR_TF).set(tf);
  integrator.input(INTEGRATOR_X0).set(x0);
  integrator.input(INTEGRATOR_P).set(u_init);

  // Integrate
  integrator.evaluate();

  cout << "before " << integrator.output().data() << endl;
  vector<double> res0 = integrator.output().data();

//  u_init += 0.01;
  x0[3] += 0.01;

  integrator.input(INTEGRATOR_X0).set(x0);
  integrator.input(INTEGRATOR_P).set(u_init);

  // Evaluate perturbed
  integrator.evaluate();
  
  // Print statistics
//   integrator.printStats();

  cout << "after " << integrator.output().data() << endl;
  
  vector<double> fd = integrator.output().data();
  for(int i=0; i<fd.size(); ++i){
    fd[i] -= res0[i];
    fd[i] /= 0.01;
  }
  
  cout << "fd    " << fd << endl;
  vector<double> x0_seed(x0.size(),0);
  double u_seed = 0;
  
  //u_seed = 1;
  x0_seed[3] = 1;

  integrator.input(INTEGRATOR_P).set(u_init);
  integrator.input(INTEGRATOR_T0).setF(0.0);
  integrator.input(INTEGRATOR_TF).setF(0.0);
  integrator.input(INTEGRATOR_X0).setF(x0_seed);
  integrator.input(INTEGRATOR_P).setF(u_seed);

  integrator.evaluate(1,0);
  vector<double> &oseed = integrator.output().dataF();
  cout << "osens " << oseed << endl;

  return 0;

  fill(oseed.begin(),oseed.end(),0);
  oseed[1] = 1;

  // Pass output seeds
  integrator.evaluate(0,1);
    
  cout << integrator.input(INTEGRATOR_T0).dataF() << endl;
  cout << integrator.input(INTEGRATOR_TF).dataF()<< endl;
  cout << integrator.input(INTEGRATOR_X0).dataF()<< endl;
  cout << integrator.input(INTEGRATOR_P).dataF() << endl;

  
  return 0;

}
