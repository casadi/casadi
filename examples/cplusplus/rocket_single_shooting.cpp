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

#include <iostream>
#include <symbolic/casadi.hpp>
#include <nonlinear_programming/scpgen.hpp>
#include <nonlinear_programming/nlp_qp_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include "interfaces/sundials/cvodes_integrator.hpp"
#include "interfaces/sundials/idas_integrator.hpp"
#include "integration/rk_integrator.hpp"

using namespace CasADi;
using namespace std;

bool sundials_integrator = true;
bool explicit_integrator = false;
bool lifted_newton = false;

int main(){
  
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
  SX t("t");

  // Differential states
  SX s("s"), v("v"), m("m");
  vector<SX> x(3); 
  x[0] = s;
  x[1] = v;
  x[2] = m;

  // Control
  SX u("u");

  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate
  
  // Differential equation
  vector<SX> rhs(3);
  rhs[0] = v;
  rhs[1] = (u-alpha*v*v)/m;
  rhs[2] = -beta*u*u;

  // Initial conditions
  vector<double> x0(3);
  x0[0] = 0;
  x0[1] = 0;
  x0[2] = 1;

  // DAE residual function
  SXFunction daefcn(daeIn<SXMatrix>("x",x, "p",u, "t",t),daeOut<SXMatrix>("ode",rhs));
  daefcn.setOption("name","DAE residual");

  // Integrator
  Integrator integrator;
  if(sundials_integrator){
    if(explicit_integrator){
      // Explicit integrator (CVODES)
      integrator = CVodesIntegrator(daefcn);
      // integrator.setOption("exact_jacobian",true);
      // integrator.setOption("linear_multistep_method","bdf"); // adams or bdf
      // integrator.setOption("nonlinear_solver_iteration","newton"); // newton or functional
    } else {
      // Implicit integrator (IDAS)
      integrator = IdasIntegrator(daefcn);
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
    integrator = RKIntegrator(daefcn);
    integrator.setOption("expand_f",true);
    integrator.setOption("interpolation_order",1);
    integrator.setOption("number_of_finite_elements",1000);
  }

  integrator.setOption("t0",t0);
  integrator.setOption("tf",tf);
  integrator.init();

  // control for all segments
  MX U("U",nu); 

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
  NLPSolver solver;
  if(lifted_newton){
    solver = SCPgen(nlp);

    solver.setOption("verbose",true);
    solver.setOption("regularize",false);
    solver.setOption("maxiter_ls",1);
    solver.setOption("maxiter",100);
    
    // Use IPOPT as QP solver
    solver.setOption("qp_solver",NLPQPSolver::creator);
    Dictionary qp_solver_options;
    qp_solver_options["nlp_solver"] = IpoptSolver::creator;
    Dictionary ipopt_options;
    ipopt_options["tol"] = 1e-12;
    ipopt_options["print_level"] = 0;
    ipopt_options["print_time"] = false;
    qp_solver_options["nlp_solver_options"] = ipopt_options;
    solver.setOption("qp_solver_options",qp_solver_options);
  } else {
    solver = IpoptSolver(nlp);
    
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
  solver.solve();

  // Get the solution
  solver.getOutput(Usol,"x");
  cout << "optimal solution: " << Usol << endl;

  return 0;
}
