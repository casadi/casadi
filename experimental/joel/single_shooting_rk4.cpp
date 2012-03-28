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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSefcn.  See the GNU
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
#include <casadi/casadi.hpp>
#include <nonlinear_programming/sqp_method.hpp>
#include <interfaces/qpoases/qpoases_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <nonlinear_programming/lifted_sqp.hpp>

using namespace CasADi;
using namespace std;

int main(){
  cout << "program started" << endl;

  // Dimensions
  int nk = 15;  // Number of control segments
  int nj = 100; // Number of integration steps per control segment

  // Control
  SXMatrix u = ssym("u",nk); // control

  // Number of states
  int nx = 3;
  
  // Intermediate variables with initial values and bounds
  SXMatrix v, v_def;
  DMatrix v_init, v_min, v_max;
  
  // Initial values and bounds for the state at the different stages
  DMatrix x_k_init =  DMatrix::zeros(nx);
  DMatrix x_k_min  = -DMatrix::inf(nx); 
  DMatrix x_k_max  =  DMatrix::inf(nx);
  
  // Initial conditions
  DMatrix x_0 = DMatrix::zeros(nx);
  x_0[0] = 0; // x
  x_0[1] = 1; // y
  x_0[2] = 0; // lterm
  
  double tf = 10;
  SX dt = tf/(nj*nk); // time step

  // For all the shooting intervals
  SXMatrix x_k = x_0;
  SXMatrix ode_rhs(x_k.sparsity(),0);
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
  SXMatrix f = x_k[2] + (tf/nk)*inner_prod(u,u);
  
  // Terminal constraints
  SXMatrix g;
  g.append(x_k[0]);
  g.append(x_k[1]);

  // Bounds on g
  DMatrix g_min = DMatrix::zeros(2);
  DMatrix g_max = DMatrix::zeros(2);

  // Bounds on u and initial condition
  DMatrix u_min  = -0.75*DMatrix::ones(nk);
  DMatrix u_max  =  1.00*DMatrix::ones(nk);
  DMatrix u_init =       DMatrix::zeros(nk);
  DMatrix xv_min = vertcat(u_min,v_min);
  DMatrix xv_max = vertcat(u_max,v_max);
  DMatrix xv_init = vertcat(u_init,v_init);
  DMatrix gv_min = vertcat(DMatrix::zeros(v.size()),g_min);
  DMatrix gv_max = vertcat(DMatrix::zeros(v.size()),g_max);
  
  // Formulate the full-space NLP
  SXFunction ffcn(vertcat(u,v),f);
  SXFunction gfcn(vertcat(u,v),vertcat(v_def-v,g));

  // Solve using IPOPT
  IpoptSolver ipopt_solver(ffcn,gfcn);
  
  // Set options
  ipopt_solver.setOption("generate_hessian",true);
  ipopt_solver.setOption("tol",1e-15);
  
  // initialize the solver
  ipopt_solver.init();

  // Initial guess and bounds
  ipopt_solver.setInput(xv_min,NLP_LBX);
  ipopt_solver.setInput(xv_max,NLP_UBX);
  ipopt_solver.setInput(xv_init,NLP_X_INIT);
  ipopt_solver.setInput(gv_min,NLP_LBG);
  ipopt_solver.setInput(gv_max,NLP_UBG);

  // Solve the problem
  ipopt_solver.solve();
  
  // Print the optimal solution
  cout << "I: optimal cost:    " << ipopt_solver.output(NLP_COST).toScalar() << endl;
  cout << "I: optimal control: " << ipopt_solver.output(NLP_X_OPT) << endl;
  cout << "I: multipliers (u): " << ipopt_solver.output(NLP_LAMBDA_X) << endl;
  cout << "I: multipliers (gb): " << ipopt_solver.output(NLP_LAMBDA_G) << endl;
  
  // Solve using lifted SQP solver
  LiftedSQP lifted_sqp(ffcn,gfcn);
  
  // Set options
  lifted_sqp.setOption("qp_solver",Interfaces::QPOasesSolver::creator);
  Dictionary qp_solver_options;
  qp_solver_options["printLevel"] = "none";
  lifted_sqp.setOption("qp_solver_options",qp_solver_options);
  lifted_sqp.setOption("num_lifted",v.size());
  
  // initialize the solver
  lifted_sqp.init();

  lifted_sqp.setInput(xv_min,NLP_LBX);
  lifted_sqp.setInput(xv_max,NLP_UBX);
  lifted_sqp.setInput(xv_init,NLP_X_INIT);
  lifted_sqp.setInput(gv_min,NLP_LBG);
  lifted_sqp.setInput(gv_max,NLP_UBG);

  // Solve the problem
  lifted_sqp.solve();
  
  // Print the optimal solution
  cout << "L: optimal cost:    " << lifted_sqp.output(NLP_COST).toScalar() << endl;
  cout << "L: optimal control: " << lifted_sqp.output(NLP_X_OPT) << endl;
  cout << "L: multipliers (u): " << lifted_sqp.output(NLP_LAMBDA_X) << endl;
  cout << "L: multipliers (gb): " << lifted_sqp.output(NLP_LAMBDA_G) << endl;
  
  // Solve using old SQP method
  SQPMethod sqp_method(ffcn,gfcn);
  
  // Set options
  sqp_method.setOption("qp_solver",Interfaces::QPOasesSolver::creator);
  sqp_method.setOption("qp_solver_options",qp_solver_options);
  sqp_method.setOption("generate_hessian",true);
  
  // initialize the solver
  sqp_method.init();

  sqp_method.setInput(xv_min,NLP_LBX);
  sqp_method.setInput(xv_max,NLP_UBX);
  sqp_method.setInput(xv_init,NLP_X_INIT);
  sqp_method.setInput(gv_min,NLP_LBG);
  sqp_method.setInput(gv_max,NLP_UBG);

  // Solve the problem
  sqp_method.solve();
  
  // Print the optimal solution
  cout << "S: optimal cost:    " << sqp_method.output(NLP_COST).toScalar() << endl;
  cout << "S: optimal control: " << sqp_method.output(NLP_X_OPT) << endl;
  cout << "S: multipliers (u): " << sqp_method.output(NLP_LAMBDA_X) << endl;
  cout << "S: multipliers (gb): " << sqp_method.output(NLP_LAMBDA_G) << endl;
  
  
  return 0;
}




