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

#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/idas_integrator.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>
#include <interfaces/sundials/kinsol_solver.hpp>
#include <interfaces/csparse/csparse.hpp>

#include <optimal_control/symbolic_ocp.hpp>
#include <optimal_control/direct_multiple_shooting.hpp>

using namespace CasADi;
using namespace std;

int main(){

  // Allocate an OCP object
  SymbolicOCP ocp;

  // Load the XML file
  ocp.parseFMI("../examples/xml_files/cstr.xml");

  // Transform into an explicit ODE
  ocp.makeExplicit();

  // Scale the variables
  ocp.scaleVariables();
  
  // Eliminate the independent parameters
  ocp.eliminateIndependentParameters();

  // Print the ocp to screen
  ocp.print();
  
  // Correct the inital guess and bounds on variables
  ocp.setStart("u",280);
  ocp.setMin("u",230);
  ocp.setMax("u",370);

  // Correct bound on state
  ocp.setMax("cstr.T",350);
      
  // Initial guess and bounds for the state
  vector<double> x0 = ocp.start(ocp.x,true);
  vector<double> xmin = ocp.min(ocp.x,true);
  vector<double> xmax = ocp.max(ocp.x,true);
  
  // Initial guess and bounds for the control
  vector<double> u0 = ocp.start(ocp.u,true);
  vector<double> umin = ocp.min(ocp.u,true);
  vector<double> umax = ocp.max(ocp.u,true);
  
  // Number of shooting nodes
  int num_nodes = 100;

  // Set integrator options
  Dictionary integrator_options;
  integrator_options["abstol"]=1e-8;
  integrator_options["reltol"]=1e-8;
  integrator_options["tf"]=ocp.tf/num_nodes;

  // Mayer objective function
  SXFunction mterm(ocp.x,ocp("cost"));
  
  // DAE residual function
  SXFunction dae(daeIn("x",ocp.x,"p",ocp.u,"t",ocp.t),daeOut("ode",ocp.ode));

  // Create a multiple shooting discretization
  DirectMultipleShooting ocp_solver;
  ocp_solver = DirectMultipleShooting(dae,mterm);
  ocp_solver.setOption("integrator",IdasIntegrator::creator);
  ocp_solver.setOption("integrator_options",integrator_options);
  ocp_solver.setOption("number_of_grid_points",num_nodes);
  ocp_solver.setOption("final_time",ocp.tf);
  ocp_solver.setOption("parallelization","openmp");
//  ocp_solver.setOption("parallelization","expand");
//  ocp_solver.setOption("parallelization","serial");

  // NLP solver
  ocp_solver.setOption("nlp_solver",IpoptSolver::creator);
  Dictionary nlp_solver_dict;
  nlp_solver_dict["tol"] = 1e-5;
  nlp_solver_dict["hessian_approximation"] = "limited-memory"; // For BFGS
  nlp_solver_dict["max_iter"] = 100;
  nlp_solver_dict["linear_solver"] = "ma57";
  //  nlp_solver_dict["derivative_test"] = "first-order";
  //  nlp_solver_dict["verbose"] = true;
  ocp_solver.setOption("nlp_solver_options",nlp_solver_dict);
  ocp_solver.init();

  // Initial condition
  for(int i=0; i<ocp.x.size(); ++i){
    ocp_solver.input("x_init")(i,0) = ocp_solver.input("lbx")(i,0) = ocp_solver.input("ubx")(i,0) = x0[i];
  }

  // State bounds
  for(int k=1; k<=num_nodes; ++k){
    for(int i=0; i<ocp.x.size(); ++i){
      ocp_solver.input("x_init")(i,k) = x0[i];
      ocp_solver.input("lbx")(i,k) = xmin[i];
      ocp_solver.input("ubx")(i,k) = xmax[i];
    }
  }

  // Control bounds
  for(int k=0; k<num_nodes; ++k){
    for(int i=0; i<ocp.u.size(); ++i){
      ocp_solver.input("u_init")(i,k) = u0[i];
      ocp_solver.input("lbu")(i,k) = umin[i];
      ocp_solver.input("ubu")(i,k) = umax[i];
    }
  }

  // Solve the problem
  ocp_solver.solve();
  
    
  
  return 0;
}
