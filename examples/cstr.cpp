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

#include <optimal_control/fmi_parser.hpp>
#include <optimal_control/fmi_parser_internal.hpp>
#include <optimal_control/ocp_tools.hpp>
#include <optimal_control/variable_tools.hpp>
#include <optimal_control/multiple_shooting.hpp>

using namespace CasADi;
using namespace CasADi::Interfaces;
using namespace CasADi::Sundials;
using namespace CasADi::OptimalControl;
using namespace std;

int main(){

  // Allocate a parser and load the xml
  FMIParser parser("../examples/xml_files/cstr.xml");

  // Obtain the symbolic representation of the OCP
  parser.parse();

  // Print the ocp to screen
  parser.print();

  // Scale the variables
  parser->scaleVariables();

  // Eliminate the dependent variables
  parser->eliminateDependent();

  // Scale the equations
  parser->scaleEquations();
  
  // Correct the inital guess and bounds on variables
  parser->u_[0].setStart(280);
  parser->u_[0].setMin(230);
  parser->u_[0].setMax(370);
  
  // Correct bound on state
  parser->x_[1].setMax(350);
  
  // Variables
  SX t = parser->t_;
  Matrix<SX> x = var(parser->x_);
  Matrix<SX> xdot = der(parser->x_);
  casadi_assert(parser->z_.empty());
  Matrix<SX> p = var(parser->p_);
  Matrix<SX> u = var(parser->u_);

  // Initial guess and bounds for the state
  vector<double> x0 = getStart(parser->x_,true);
  vector<double> xmin = getMin(parser->x_,true);
  vector<double> xmax = getMax(parser->x_,true);
  
  // Initial guess and bounds for the control
  vector<double> u0 = getStart(parser->u_,true);
  vector<double> umin = getMin(parser->u_,true);
  vector<double> umax = getMax(parser->u_,true);
  
  // Integrator instance
  Integrator integrator;

  // Create an implicit function residual
  vector<Matrix<SX> > impres_in(DAE_NUM_IN+1);
  impres_in[0] = xdot;
  impres_in[1+DAE_T] = t;
  impres_in[1+DAE_Y] = x;
  impres_in[1+DAE_P] = u;
  SXFunction impres(impres_in,parser->implicit_fcn_);
  
  // Create an implicit function (KINSOL)
  KinsolSolver ode(impres);
  ode.setLinearSolver(CSparse(CRSSparsity()));
  ode.setOption("linear_solver","user_defined");
  ode.init();
  
  // DAE residual
  vector<Matrix<SX> > dae_in(DAE_NUM_IN);
  dae_in[DAE_T] = t;
  dae_in[DAE_Y] = x;
  dae_in[DAE_YDOT] = xdot;
  dae_in[DAE_P] = u;
  SXFunction dae(dae_in,parser->implicit_fcn_);

  bool use_kinsol = false;
  if(use_kinsol){
    // Create an ODE integrator (CVodes)
    integrator = CVodesIntegrator(ode);
    
  } else {
    // Create DAE integrator (IDAS)
    integrator = IdasIntegrator(dae);
    
  }

  // Number of shooting nodes
  int num_nodes = 100;

  // Set integrator options
  integrator.setOption("number_of_fwd_dir",1);
  integrator.setOption("number_of_adj_dir",0);
  integrator.setOption("exact_jacobian",true);
  integrator.setOption("fsens_err_con",true);
  integrator.setOption("quad_err_con",true);
  integrator.setOption("abstol",1e-8);
  integrator.setOption("reltol",1e-8);
  integrator.setOption("store_jacobians",true);
  integrator.setOption("tf",parser->tf/num_nodes);
  integrator.init();

  // Mayer objective function
  Matrix<SX> xf = symbolic("xf",x.size(),1);
  SXFunction mterm(xf, xf[0]);
  mterm.setOption("store_jacobians",true);
  
  // Create a multiple shooting discretization
  MultipleShooting ms(integrator,mterm);
  ms.setOption("number_of_grid_points",num_nodes);
  ms.setOption("final_time",parser->tf);
  ms.setOption("parallelization","openmp");
//  ms.setOption("parallelization","expand");

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

  // Initial condition
  for(int i=0; i<x.size(); ++i){
    ms.input(OCP_X_INIT)(i,0) = ms.input(OCP_LBX)(i,0) = ms.input(OCP_UBX)(i,0) = x0[i];
  }

  // State bounds
  for(int k=1; k<=num_nodes; ++k){
    for(int i=0; i<x.size(); ++i){
      ms.input(OCP_X_INIT)(i,k) = x0[i];
      ms.input(OCP_LBX)(i,k) = xmin[i];
      ms.input(OCP_UBX)(i,k) = xmax[i];
    }
  }

  // Control bounds
  for(int k=0; k<num_nodes; ++k){
    for(int i=0; i<u.size(); ++i){
      ms.input(OCP_U_INIT)(i,k) = u0[i];
      ms.input(OCP_LBU)(i,k) = umin[i];
      ms.input(OCP_UBU)(i,k) = umax[i];
    }
  }

  // Solve the problem
  ms.solve();
  
    
  
  return 0;
}
