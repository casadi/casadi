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

#include <casadi/stl_vector_tools.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/idas_integrator.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>
#include <interfaces/sundials/kinsol_solver.hpp>
#include <interfaces/csparse/csparse.hpp>

#include <casadi/fx/fx_tools.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>
#include <casadi/fx/jacobian.hpp>

#include <optimal_control/symbolic_ocp.hpp>
#include <optimal_control/ocp_tools.hpp>
#include <optimal_control/variable_tools.hpp>
#include <optimal_control/multiple_shooting.hpp>

using namespace CasADi;
using namespace CasADi::Interfaces;
using namespace CasADi::Sundials;
using namespace CasADi::OptimalControl;
using namespace std;

int main(){
  bool use_kinsol = false;

  // Allocate an OCP object
  SymbolicOCP ocp;

  // Load the XML file
  Dictionary parse_options;
  parse_options["scale_variables"] = true;
  parse_options["eliminate_dependent"] = true;
  parse_options["scale_equations"] = false;
  ocp.parseFMI("../examples/xml_files/cstr.xml",parse_options);
  
  // To explicit form
//   ocp.makeExplicit();
  
  // Print the ocp to screen
  ocp.print();
  
  // Correct the inital guess and bounds on variables
  ocp.variable("u").setStart(280);
  ocp.variable("u").setMin(230);
  ocp.variable("u").setMax(370);

  // Correct bound on state
  ocp.variable("cstr.T").setMax(350);
  
  // Variables
  SXMatrix t = ocp.t;
  SXMatrix x = var(ocp.x);
  SXMatrix xdot = der(ocp.x);
  SXMatrix p = var(ocp.p);
  SXMatrix u = var(ocp.u);
    
  // Initial guess and bounds for the state
  vector<double> x0 = getStart(ocp.x,true);
  vector<double> xmin = getMin(ocp.x,true);
  vector<double> xmax = getMax(ocp.x,true);
  
  // Initial guess and bounds for the control
  vector<double> u0 = getStart(ocp.u,true);
  vector<double> umin = getMin(ocp.u,true);
  vector<double> umax = getMax(ocp.u,true);
  
  // Number of shooting nodes
  int num_nodes = 100;

  // Set integrator options
  Dictionary integrator_options;
  integrator_options["number_of_fwd_dir"]=1;
  integrator_options["number_of_adj_dir"]=0;
  integrator_options["exact_jacobian"]=true;
  integrator_options["fsens_err_con"]=true;
  integrator_options["quad_err_con"]=true;
  integrator_options["abstol"]=1e-8;
  integrator_options["reltol"]=1e-8;
  integrator_options["store_jacobians"]=true;
  integrator_options["tf"]=ocp.tf/num_nodes;

  // Mayer objective function
  SXMatrix xf = ssym("xf",x.size(),1);
  SXFunction mterm(xf, xf[0]);
  mterm.setOption("store_jacobians",true);
  
  // Create a multiple shooting discretization
  MultipleShooting ocp_solver;
  if(use_kinsol){
    // Create an implicit function residual
    vector<SXMatrix > impres_in(DAE_NUM_IN+1);
    impres_in[0] = xdot;
    impres_in[1+DAE_T] = t;
    impres_in[1+DAE_X] = x;
    impres_in[1+DAE_P] = u;
    SXFunction impres(impres_in,ocp.dae);

    // Create an implicit function (KINSOL)
    KinsolSolver ode(impres);
    ode.setOption("linear_solver_creator",CSparse::creator);
    ode.init();
    
    // Create OCP solver
    ocp_solver = MultipleShooting(ode,mterm);
    ocp_solver.setOption("integrator",CVodesIntegrator::creator);
  } else {
    // DAE residual function
    SXFunction dae(daeIn(x,SXMatrix(),u,t,xdot),daeOut(ocp.dae));
    
    ocp_solver = MultipleShooting(dae,mterm);
    ocp_solver.setOption("integrator",IdasIntegrator::creator);
  }
  ocp_solver.setOption("integrator_options",integrator_options);
  ocp_solver.setOption("number_of_grid_points",num_nodes);
  ocp_solver.setOption("final_time",ocp.tf);
  ocp_solver.setOption("parallelization","openmp");
//  ocp_solver.setOption("parallelization","expand");

  // NLP solver
  ocp_solver.setOption("nlp_solver",IpoptSolver::creator);
  Dictionary nlp_solver_dict;
  nlp_solver_dict["tol"] = 1e-5;
  nlp_solver_dict["hessian_approximation"] = "limited-memory";
  nlp_solver_dict["max_iter"] = 100;
  nlp_solver_dict["linear_solver"] = "ma57";
  //  nlp_solver_dict["derivative_test"] = "first-order";
  //  nlp_solver_dict["verbose"] = true;
  ocp_solver.setOption("nlp_solver_options",nlp_solver_dict);
  ocp_solver.init();

  // Initial condition
  for(int i=0; i<x.size(); ++i){
    ocp_solver.input(OCP_X_INIT)(i,0) = ocp_solver.input(OCP_LBX)(i,0) = ocp_solver.input(OCP_UBX)(i,0) = x0[i];
  }

  // State bounds
  for(int k=1; k<=num_nodes; ++k){
    for(int i=0; i<x.size(); ++i){
      ocp_solver.input(OCP_X_INIT)(i,k) = x0[i];
      ocp_solver.input(OCP_LBX)(i,k) = xmin[i];
      ocp_solver.input(OCP_UBX)(i,k) = xmax[i];
    }
  }

  // Control bounds
  for(int k=0; k<num_nodes; ++k){
    for(int i=0; i<u.size(); ++i){
      ocp_solver.input(OCP_U_INIT)(i,k) = u0[i];
      ocp_solver.input(OCP_LBU)(i,k) = umin[i];
      ocp_solver.input(OCP_UBU)(i,k) = umax[i];
    }
  }

  // Solve the problem
  ocp_solver.solve();
  
    
  
  return 0;
}
