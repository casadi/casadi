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

  // Correct the inital guess and bounds on variables
  ocp.setStart("u",280);
  ocp.setMin("u",230);
  ocp.setMax("u",370);

  // Correct bound on state
  ocp.setMax("cstr.T",350);

  // Transform into an explicit ODE
  ocp.makeExplicit();

  // Scale the variables
  ocp.scaleVariables();
  
  // Eliminate the independent parameters
  ocp.eliminateIndependentParameters();

  // Print the ocp to screen
  ocp.print();
        
  // Initial guess 

  // Function to evaluate the initial guess and bounds for all time points
  vector<SX> bfcn_out;
  bfcn_out.push_back(ocp.start(ocp.x));
  bfcn_out.push_back(ocp.min(ocp.x));
  bfcn_out.push_back(ocp.max(ocp.x));
  bfcn_out.push_back(ocp.start(ocp.u));
  bfcn_out.push_back(ocp.min(ocp.u));
  bfcn_out.push_back(ocp.max(ocp.u));
  SXFunction bfcn(ocp.t,bfcn_out);
  bfcn.init();
  
  // Number of shooting nodes
  int num_nodes = 100;

  // Set integrator options
  Dictionary integrator_options;
  integrator_options["abstol"]=1e-8;
  integrator_options["reltol"]=1e-8;
  integrator_options["tf"]=ocp.tf/num_nodes;

  // Mayer objective function
  SXFunction mterm(ocp.x,ocp.mterm);
  
  // DAE residual function
  SXFunction dae(daeIn("x",ocp.x,"p",ocp.u,"t",ocp.t),daeOut("ode",ocp.ode(ocp.x)));

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

  // Pass the bounds and initial guess
  for(int i=0; i<num_nodes+1; ++i){
    // Evaluate the function
    bfcn.setInput((i*ocp.tf)/num_nodes);
    bfcn.evaluate();

    if(i==0){
      ocp_solver.input("x_init")(Slice(),i) = ocp_solver.input("lbx")(Slice(),i) = ocp_solver.input("ubx")(Slice(),i) = bfcn.output(0);
    } else {
      ocp_solver.input("x_init")(Slice(),i) = bfcn.output(0);
      ocp_solver.input("lbx")(Slice(),i) = bfcn.output(1);
      ocp_solver.input("ubx")(Slice(),i) = bfcn.output(2);
    }

    if(i<num_nodes){
      ocp_solver.input("u_init")(Slice(),i) = bfcn.output(3);
      ocp_solver.input("lbu")(Slice(),i) = bfcn.output(4);
      ocp_solver.input("ubu")(Slice(),i) = bfcn.output(5);
    }
  }

  // Solve the problem
  ocp_solver.solve();
  
    
  
  return 0;
}
