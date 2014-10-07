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
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std;

// Declare solvers to be loaded manually
extern "C" void casadi_load_integrator_cvodes();
extern "C" void casadi_load_integrator_idas();
extern "C" void casadi_load_integrator_rk();
extern "C" void casadi_load_nlpsolver_ipopt();
extern "C" void casadi_load_nlpsolver_scpgen();

bool sundials_integrator = true;
bool explicit_integrator = false;
bool lifted_newton = false;

int main(){
  // Load integrators manually
  casadi_load_integrator_cvodes();
  casadi_load_integrator_idas();
  casadi_load_integrator_rk();
  casadi_load_nlpsolver_ipopt();
  casadi_load_nlpsolver_scpgen();
  
  // Time length
  double T = 10.0;

  // Shooting length
  int nu = 20; // Number of control segments

  // Time horizon for integrator
  double t0 = 0;
  double tf = T/nu;

  // Initial position
  vector<double> X0(3);
  X0[0] = 0; // initial position
  X0[1] = 0; // initial speed
  X0[2] = 1; // initial mass

  // Time 
  SX t = SX::sym("t");

  // Differential states
  SX s = SX::sym("s"), v = SX::sym("v"), m = SX::sym("m");
  SX x = SX::zeros(3);
  x[0] = s;
  x[1] = v;
  x[2] = m;

  // Control
  SX u = SX::sym("u");

  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate
  
  // Differential equation
  SX rhs = SX::zeros(3);
  rhs[0] = v;
  rhs[1] = (u-alpha*v*v)/m;
  rhs[2] = -beta*u*u;

  // Initial conditions
  vector<double> x0(3);
  x0[0] = 0;
  x0[1] = 0;
  x0[2] = 1;

  // DAE residual function
  SXFunction daefcn(daeIn("x",x, "p",u, "t",t),daeOut("ode",rhs));
  daefcn.setOption("name","DAE residual");

  // Integrator
  Integrator integrator;
  if(sundials_integrator){
    if(explicit_integrator){
      // Explicit integrator (CVODES)
      integrator = Integrator("cvodes", daefcn);
      // integrator.setOption("exact_jacobian",true);
      // integrator.setOption("linear_multistep_method","bdf"); // adams or bdf
      // integrator.setOption("nonlinear_solver_iteration","newton"); // newton or functional
    } else {
      // Implicit integrator (IDAS)
      integrator = Integrator("idas", daefcn);
      integrator.setOption("calc_ic",false);
    }
    integrator.setOption("fsens_err_con",true);
    integrator.setOption("quad_err_con",true);
    integrator.setOption("abstol",1e-6);
    integrator.setOption("reltol",1e-6);
    integrator.setOption("stop_at_end",false);
    //  integrator.setOption("fsens_all_at_once",false);
    integrator.setOption("steps_per_checkpoint",100); // BUG: Too low number causes segfaults
  } else {
    // An explicit Euler integrator
    integrator = Integrator("rk", daefcn);
    integrator.setOption("expand_f",true);
    integrator.setOption("interpolation_order",1);
    integrator.setOption("number_of_finite_elements",1000);
  }

  integrator.setOption("t0",t0);
  integrator.setOption("tf",tf);
  integrator.init();

  // control for all segments
  MX U = MX::sym("U",nu); 

  // Integrate over all intervals
  MX X=X0;
  for(int k=0; k<nu; ++k){
    // Assemble the input
    vector<MX> input(INTEGRATOR_NUM_IN);
    input[INTEGRATOR_X0] = X;
    input[INTEGRATOR_P] = U[k];

    // Integrate
    X = integrator.call(input).at(0);

    // Lift X
    if(lifted_newton){
      X.lift(X);
    }
  }

  // Objective function
  MX F = inner_prod(U,U);

  // Terminal constraints
  MX G = vertcat(X[0],X[1]);
  
  // Create the NLP
  MXFunction nlp(nlpIn("x",U),nlpOut("f",F,"g",G));

  // Allocate an NLP solver
  NlpSolver solver;
  if(lifted_newton){
    solver = NlpSolver("scpgen", nlp);

    solver.setOption("verbose",true);
    solver.setOption("regularize",false);
    solver.setOption("max_iter_ls",1);
    solver.setOption("max_iter",100);
    
    // Use IPOPT as QP solver
    solver.setOption("qp_solver","nlp");
    Dictionary qp_solver_options;
    qp_solver_options["nlp_solver"] = "ipopt";
    Dictionary ipopt_options;
    ipopt_options["tol"] = 1e-12;
    ipopt_options["print_level"] = 0;
    ipopt_options["print_time"] = false;
    qp_solver_options["nlp_solver_options"] = ipopt_options;
    solver.setOption("qp_solver_options",qp_solver_options);
  } else {
    solver = NlpSolver("ipopt", nlp);
    
    // Set options
    solver.setOption("tol",1e-10);
    solver.setOption("hessian_approximation","limited-memory");
  }

  // initialize the solver
  solver.init();

  // Bounds on u and initial condition
  vector<double> Umin(nu), Umax(nu), Usol(nu);
  for(int i=0; i<nu; ++i){
    Umin[i] = -10;
    Umax[i] =  10;
    Usol[i] = 0.4;
  }
  solver.setInput(Umin,"lbx");
  solver.setInput(Umax,"ubx");
  solver.setInput(Usol,"x0");

  // Bounds on g
  vector<double> Gmin(2), Gmax(2);
  Gmin[0] = Gmax[0] = 10;
  Gmin[1] = Gmax[1] =  0;
  solver.setInput(Gmin,"lbg");
  solver.setInput(Gmax,"ubg");

  // Solve the problem
  solver.evaluate();

  // Get the solution
  solver.getOutput(Usol,"x");
  cout << "optimal solution: " << Usol << endl;

  return 0;
}
