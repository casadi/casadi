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

/** \brief Writing a multiple shooting code from scratch 

  This example demonstrates how to write a simple, yet powerful multiple shooting code from 
  scratch using CasADi. For clarity, the code below uses a simple formulation with only states 
  and controls, and no path constrants. It relies on CasADi machinery to keep track of sparsity,
  formulate ODE sensitivity equations and build up the Jacobian of the NLP constraint function.
    
  By extending the code below, it should be possible for a user to solve ever more complex 
  problems. For example, one can easily make a multi-stage formulation by simply allocating 
  another integrator instance and use the two integrators instances for different shooting nodes.
  By replacing the explicit CVodes integrator with the fully-implicit IDAS integrator, one can
  solve optimal control problems in differential-algebraic equations.
  
  \author Joel Andersson
  \date 2011-2012
*/

// CasADi core
#include <symbolic/casadi.hpp>

// Interfaces
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>

using namespace CasADi;
using namespace std;

// Infinity
double inf = numeric_limits<double>::infinity();

int main(){

  // Declare variables
  SXMatrix u = ssym("u"); // control
  SXMatrix r = ssym("r"), s = ssym("s"); // states
  SXMatrix x = vertcat(r,s);

  // Number of differential states
  int nx = x.size1();
  
  // Number of controls
  int nu = u.size1();

  // Bounds and initial guess for the control
  double u_min[] =  { -0.75 };
  double u_max[]  = {  1.0  };
  double u_init[] = {  0.0  };

  // Bounds and initial guess for the state
  double x0_min[] = {   0,    1 };
  double x0_max[] = {   0,    1 };
  double x_min[]  = {-inf, -inf };
  double x_max[]  = { inf,  inf };
  double xf_min[] = {   0,    0 };
  double xf_max[] = {   0,    0 };
  double x_init[] = {   0,    0 };

  // Final time
  double tf = 20.0;
  
  // Number of shooting nodes
  int ns = 50;

  // ODE right hand side and quadrature
  SXMatrix ode = vertcat((1 - s*s)*r - s + u, r);
  SXMatrix quad = r*r + s*s + u*u;
  SXFunction rhs(daeIn("x",x,"p",u),daeOut("ode",ode,"quad",quad));

  // Create an integrator (CVodes)
  CVodesIntegrator integrator(rhs);
  integrator.setOption("t0",0);
  integrator.setOption("tf",tf/ns);
  integrator.init();
  
  // Total number of NLP variables
  int NV = nx*(ns+1) + nu*ns;
  
  // Declare variable vector for the NLP
  MX V = msym("V",NV);

  // NLP variable bounds and initial guess
  vector<double> v_min,v_max,v_init;
  
  // Offset in V
  int offset=0; 

  // State at each shooting node and control for each shooting interval
  vector<MX> X, U;
  for(int k=0; k<ns; ++k){
    // Local state
    X.push_back( V[Slice(offset,offset+nx)] );
    if(k==0){
      v_min.insert(v_min.end(),x0_min,x0_min+nx);
      v_max.insert(v_max.end(),x0_max,x0_max+nx);
    } else {
      v_min.insert(v_min.end(),x_min,x_min+nx);
      v_max.insert(v_max.end(),x_max,x_max+nx);
    }
    v_init.insert(v_init.end(),x_init,x_init+nx);    
    offset += nx;
    
    // Local control
    U.push_back( V[Slice(offset,offset+nu)] );
    v_min.insert(v_min.end(),u_min,u_min+nu);
    v_max.insert(v_max.end(),u_max,u_max+nu);
    v_init.insert(v_init.end(),u_init,u_init+nu);    
    offset += nu;
  }
  
  // State at end
  X.push_back(V[Slice(offset,offset+nx)]);
  v_min.insert(v_min.end(),xf_min,xf_min+nx);
  v_max.insert(v_max.end(),xf_max,xf_max+nx);
  v_init.insert(v_init.end(),x_init,x_init+nx);    
  offset += nx;
  
  // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
  casadi_assert(offset==NV);

  // Objective function
  MX J = 0;
  
  //Constraint function and bounds
  vector<MX> g;

  // Loop over shooting nodes
  for(int k=0; k<ns; ++k){
    // Input to the integrator
    vector<MX> int_in(INTEGRATOR_NUM_IN);
    int_in[INTEGRATOR_P] = U[k];
    int_in[INTEGRATOR_X0] = X[k];

    // Create an evaluation node
    vector<MX> I_out = integrator.call(int_in);

    // Save continuity constraints
    g.push_back( I_out[INTEGRATOR_XF] - X[k+1] );
    
    // Add objective function contribution
    J += I_out[INTEGRATOR_QF];
  }
  
  // NLP 
  MXFunction nlp(nlpIn("x",V),nlpOut("f",J,"g",vertcat(g)));
  
  // Create an NLP solver instance
  IpoptSolver nlp_solver(nlp);
  nlp_solver.setOption("tol",1e-5);
  nlp_solver.setOption("max_iter",100);
  nlp_solver.setOption("linear_solver","ma57");
  nlp_solver.init();
    
  // Initial guess and bounds on variables
  nlp_solver.setInput(v_init,"x0");
  nlp_solver.setInput(v_min,"lbx");
  nlp_solver.setInput(v_max,"ubx");
  
  // All nonlinear constraints are equality constraints
  nlp_solver.setInput(0.,"lbg");
  nlp_solver.setInput(0.,"ubg");
  
  // Solve the problem
  nlp_solver.solve();

  // Optimal solution of the NLP
  const Matrix<double>& V_opt = nlp_solver.output("x");
  
  // Get the optimal state trajectory
  vector<double> r_opt(ns+1), s_opt(ns+1);
  for(int i=0; i<=ns; ++i){
    r_opt[i] = V_opt.at(i*(nx+1));
    s_opt[i] = V_opt.at(1+i*(nx+1));
  }
  cout << "r_opt = " << endl << r_opt << endl;
  cout << "s_opt = " << endl << s_opt << endl;
  
  // Get the optimal control
  vector<double> u_opt(ns);
  for(int i=0; i<ns; ++i){
    u_opt[i] = V_opt.at(nx + i*(nx+1));
  }
  cout << "u_opt = " << endl << u_opt << endl;
  
  
  return 0;
}
