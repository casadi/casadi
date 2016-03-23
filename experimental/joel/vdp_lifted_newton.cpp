/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <core/casadi.hpp>
#include <nonlinear_programming/sqp_method.hpp>
#include <interfaces/qpoases/qpoases_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <core/std_vector_tools.hpp>
#include <nonlinear_programming/lifted_sqp.hpp>

using namespace casadi;
using namespace std;

int main(){
  cout << "program started" << endl;

  // Dimensions
  int nk = 100;  // Number of control segments
  int nj = 100; // Number of integration steps per control segment

  // Control
  SX u = ssym("u",nk); // control

  // Number of states
  int nx = 3;
  
  // Intermediate variables with initial values and bounds
  SX v, v_def;
  DM v_init, v_min, v_max;
  
  // Initial values and bounds for the state at the different stages
  DM x_k_init =  DM::zeros(nx);
  DM x_k_min  = -DM::inf(nx); 
  DM x_k_max  =  DM::inf(nx);
  
  // Initial conditions
  DM x_0 = DM::zeros(nx);
  x_0[0] = 0; // x
  x_0[1] = 1; // y
  x_0[2] = 0; // lterm
  
  double tf = 10;
  SX dt = tf/(nj*nk); // time step

  // For all the shooting intervals
  SX x_k = x_0;
  SX ode_rhs(x_k.sparsity(),0);
  for(int k=0; k<nk; ++k){
    // Get control
    SX u_k = u[k].at(0);
    
    // Integrate over the interval with Euler forward
    for(int j=0; j<nj; ++j){
      
      // ODE right hand side
      ode_rhs[0] = (1 - x_k[1]*x_k[1])*x_k[0] - x_k[1] + u_k;
      ode_rhs[1] = x_k[0];
      ode_rhs[2] = x_k[0]*x_k[0] + x_k[1]*x_k[1];
      
      // Take a step
      x_k += dt*ode_rhs;
    }
    
    // Lift x
    v_def.append(x_k);
    v_init.append(x_k_init);
    v_min.append(x_k_min);
    v_max.append(x_k_max);
    
    // Allocate intermediate variables
    stringstream ss;
    ss << "v_" << k;
    x_k = ssym(ss.str(),nx);
    v.append(x_k);
  }

  // Objective function
  SX f = x_k[2] + (tf/nk)*dot(u,u);
  
  // Terminal constraints
  SX g;
  g.append(x_k[0]);
  g.append(x_k[1]);

  // Bounds on g
  DM g_min = DM::zeros(2);
  DM g_max = DM::zeros(2);

  // Bounds on u and initial condition
  DM u_min  = -0.75*DM::ones(nk);
  DM u_max  =  1.00*DM::ones(nk);
  DM u_init =       DM::zeros(nk);
  DM xv_min = vertcat(u_min,v_min);
  DM xv_max = vertcat(u_max,v_max);
  DM xv_init = vertcat(u_init,v_init);
  DM gv_min = vertcat(DM::zeros(v.size()),g_min);
  DM gv_max = vertcat(DM::zeros(v.size()),g_max);
  
  // Formulate the full-space NLP
  SXFunction ffcn(vertcat(u,v),f);
  SXFunction gfcn(vertcat(u,v),vertcat(v_def-v,g));

  Dictionary qpsol_options;
  qpsol_options["printLevel"] = "none";
  
  // Solve using multiple NLP solvers
  enum Tests{IPOPT, LIFTED_SQP, FULLSPACE_SQP, OLD_SQP_METHOD, NUM_TESTS};
  for(int test=0; test<NUM_TESTS; ++test){
    // Get the nlp solver and NLP solver options
    NlpSolver nlp_solver;
    switch(test){
      case IPOPT:
	cout << "Testing IPOPT" << endl;
	nlp_solver = IpoptSolver(ffcn,gfcn);
	nlp_solver.setOption("generate_hessian",true);
	nlp_solver.setOption("tol",1e-10);
	break;
      case LIFTED_SQP:
	cout << "Testing lifted SQP" << endl;
	nlp_solver = LiftedSQP(ffcn,gfcn);
	nlp_solver.setOption("qpsol",QPOasesSolver::creator);
	nlp_solver.setOption("qpsol_options",qpsol_options);
	nlp_solver.setOption("num_lifted",v.size());
	nlp_solver.setOption("toldx",1e-10);
	nlp_solver.setOption("verbose",true);
	break;
      case FULLSPACE_SQP:
	cout << "Testing fullspace SQP" << endl;
	nlp_solver = LiftedSQP(ffcn,gfcn);
	nlp_solver.setOption("qpsol",QPOasesSolver::creator);
	nlp_solver.setOption("qpsol_options",qpsol_options);
	nlp_solver.setOption("num_lifted",0);
	nlp_solver.setOption("toldx",1e-10);
	nlp_solver.setOption("verbose",true);
	break;
      case OLD_SQP_METHOD:
	cout << "Testing old SQP method" << endl;
	nlp_solver = SQPMethod(ffcn,gfcn);
	nlp_solver.setOption("qpsol",QPOasesSolver::creator);
	nlp_solver.setOption("qpsol_options",qpsol_options);
	nlp_solver.setOption("generate_hessian",true);
    }
    
    // initialize the solver
    nlp_solver.init();

    // Initial guess and bounds
    nlp_solver.setInput(xv_min,"lbx");
    nlp_solver.setInput(xv_max,"ubx");
    nlp_solver.setInput(xv_init,"x0");
    nlp_solver.setInput(gv_min,"lbg");
    nlp_solver.setInput(gv_max,"ubg");

    // Solve the problem
    nlp_solver.solve();
    
    // Print the optimal solution
//     cout << "optimal cost:    " << nlp_solver.output(NLP_SOLVER_F).scalar() << endl;
//     cout << "optimal control: " << nlp_solver.output(NLP_SOLVER_X) << endl;
//     cout << "multipliers (u): " << nlp_solver.output(NLP_SOLVER_LAM_X) << endl;
//     cout << "multipliers (gb): " << nlp_solver.output(NLP_SOLVER_LAM_G) << endl;
  }
  
  return 0;
}




