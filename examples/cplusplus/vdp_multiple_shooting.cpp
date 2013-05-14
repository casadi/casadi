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

#include <symbolic/casadi.hpp>
                                   
// Solvers
#include <integration/collocation_integrator.hpp>
#include <nonlinear_programming/newton_implicit_solver.hpp>
#include <optimal_control/direct_multiple_shooting.hpp>

// 3rd party interfaces
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>
#include <interfaces/sundials/idas_integrator.hpp>
#include <interfaces/sundials/kinsol_solver.hpp>
#include <interfaces/csparse/csparse.hpp>

bool use_collocation_integrator = false;

using namespace CasADi;
using namespace std;


int main(){
  //Final time (fixed)
  double tf = 10.0;

  // Infinity
  double inf = numeric_limits<double>::infinity();

  // Declare variables (use simple, efficient DAG)
  SXMatrix x = ssym("x");
  SXMatrix y = ssym("y");
  SXMatrix u = ssym("u");
  SXMatrix L = ssym("cost");
  
  // All states
  SXMatrix states = SXMatrix::zeros(3);
  states(0) = x;
  states(1) = y;
  states(2) = L;

  //ODE right hand side
  SXMatrix f = SXMatrix::zeros(3);
  f(0) = (1 - y*y)*x - y + u;
  f(1) = x;
  f(2) = x*x + y*y + u*u;
  
  // DAE residual
  SXFunction res(daeIn("x",states, "p",u),daeOut("ode",f));
  
  Dictionary integrator_options;
  if(use_collocation_integrator){
    // integrator_options["implicit_solver"] = KinsolSolver::creator;    
    integrator_options["implicit_solver"] = NewtonImplicitSolver::creator;
    Dictionary implicit_solver_options;
    implicit_solver_options["linear_solver"] = CSparse::creator;
    integrator_options["implicit_solver_options"] = implicit_solver_options;
  } else {
    integrator_options["abstol"]=1e-8; //abs. tolerance
    integrator_options["reltol"]=1e-8; //rel. tolerance
    integrator_options["steps_per_checkpoint"]=500;
    integrator_options["stop_at_end"]=true;
  }
  
  //Numboer of shooting nodes
  int ns = 50;

  // Number of differential states
  int nx = 3;
  
  // Number of controls
  int nu = 1;
  
  // Mayer objective function
  Matrix<SX> xf = ssym("xf",nx,1);
  SXFunction mterm(xf, xf[nx-1]);

  // Create a multiple shooting discretization
  DirectMultipleShooting ms(res,mterm);
  if(use_collocation_integrator){
    ms.setOption("integrator",CollocationIntegrator::creator);
  } else {
    ms.setOption("integrator",CVodesIntegrator::creator);
    //ms.setOption("integrator",IdasIntegrator::creator);
  }
  ms.setOption("integrator_options",integrator_options);
  ms.setOption("number_of_grid_points",ns);
  ms.setOption("final_time",tf);
  ms.setOption("parallelization","openmp");
  
  // NLP solver
  ms.setOption("nlp_solver",IpoptSolver::creator);
  Dictionary nlp_solver_dict;
  nlp_solver_dict["tol"] = 1e-5;
  nlp_solver_dict["hessian_approximation"] = "limited-memory";
  nlp_solver_dict["max_iter"] = 100;
  nlp_solver_dict["linear_solver"] = "ma57";
  //  nlp_solver_dict["derivative_test"] = "first-order";
  //  nlp_solver_dict["verbose"] = true;
  ms.setOption("nlp_solver_options",nlp_solver_dict);
  
  ms.init();

  //Control bounds
  double u_min[] = {-0.75};
  double u_max[] = {1.0};
  double u_init[] = {0.0};
  
  for(int k=0; k<ns; ++k){
    copy(u_min,u_min+nu,ms.input("lbu").begin()+k*nu);
    copy(u_max,u_max+nu,ms.input("ubu").begin()+k*nu);
    copy(u_init,u_init+nu,ms.input("u_init").begin()+k*nu);
  }
  
  ms.input("lbx").setAll(-inf);
  ms.input("ubx").setAll(inf);
  ms.input("x_init").setAll(0);

  // Initial condition
  ms.input("lbx")(0,0) = ms.input("ubx")(0,0) = 0;
  ms.input("lbx")(1,0) = ms.input("ubx")(1,0) = 1;
  ms.input("lbx")(2,0) = ms.input("ubx")(2,0) = 0;

  // Final condition
  ms.input("lbx")(0,ns) = ms.input("ubx")(0,ns) = 0; 
  ms.input("lbx")(1,ns) = ms.input("ubx")(1,ns) = 0; 
  
  // Solve the problem
  ms.solve();

  cout << ms.output("x_opt") << endl;
  cout << ms.output("u_opt") << endl;
  
  return 0;
}



