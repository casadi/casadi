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

/** Writing a multiple shooting code from scratch 

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
  \date 2011
*/



// CasADi core
#include <casadi/fx/fx_tools.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/fx/sx_function.hpp>
#include <casadi/fx/mx_function.hpp>

// Interfaces
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>

using namespace CasADi;
using namespace CasADi::Sundials;
using namespace std;

// Infinity
double inf = numeric_limits<double>::infinity();

// Exact hessian?
bool exact_hessian = false;

int main(){

  // Declare variables
  SX t("t"); // time
  SX u("u"); // control
  SX r("r"), s("s"), lterm("lterm"); // states

  // Number of differential states
  const int nx = 3;
  
  // Number of controls
  const int nu = 1;

  // All states
  vector<SX> x(nx);  x[0] = r;  x[1] = s;  x[2] = lterm;

  // Bounds and initial guess for the control
  double u_min[nu] =  { -0.75 };
  double u_max[nu]  = {  1.0  };
  double u_init[nu] = {  0.0  };

  // Bounds and initial guess for the state
  double x0_min[nx] = {   0,   1,   0 };
  double x0_max[nx] = {   0,   1,   0 };
  double x_min[nx]  = {-inf,-inf,-inf };
  double x_max[nx]  = { inf, inf, inf };
  double xf_min[nx] = {   0,   0,-inf };
  double xf_max[nx] = {   0,   0, inf };
  double x_init[nx] = {   0,   0,   0 };

  // Final time
  double tf = 20.0;
  
  // Number of shooting nodes
  int ns = 50;

  // Input to the ODE/DAE functions
  vector<SXMatrix> rhs_in = daeIn<SXMatrix>(SXMatrix(),x,SXMatrix(),u,t);

  // ODE right hand side
  vector<SX> f(3);
  f[0] = (1 - s*s)*r - s + u;
  f[1] = r;
  f[2] = r*r + s*s + u*u;
  SXFunction rhs(rhs_in,daeOut<SXMatrix>(f));

  // Mayer objective function
  SXFunction mterm(x, lterm);
  mterm.init();

  //Create an integrator (CVodes)
  CVodesIntegrator integrator(rhs);
  integrator.setOption("abstol",1e-8); //abs. tolerance
  integrator.setOption("reltol",1e-8); //rel. tolerance
  integrator.setOption("steps_per_checkpoint",500);
  integrator.setOption("stop_at_end",true);
  integrator.setOption("t0",0);
  integrator.setOption("tf",tf/ns);
  integrator.setOption("number_of_fwd_dir",7);
  //integrator.setOption("verbose",true);
  integrator.init();
  
  // Total number of NLP variables
  int NV = nx*(ns+1) + nu*ns;
           
  // Declare variable vector for the NLP
  MX V("V",NV);

  // offset in the variable vector
  int el=0; 
  
  // Disretized variables for each shooting node
  vector<MX> X(ns+1), U(ns);
  for(int k=0; k<=ns; ++k){ // interior nodes
    // Local state
    X[k] = V[range(el,el+nx)];
    el += nx;
    
    // Variables below do not appear at the end point
    if(k==ns) break;
    
    // Local control
    U[k] = V[range(el,el+nu)];
    el += nu;
  }
  
  // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
  casadi_assert(el==NV);

  //Constraint function
  vector<MX> g(ns);

  // Loop over shooting nodes
  for(int k=0; k<ns; ++k){
    // Input to the integrator
    vector<MX> int_in(INTEGRATOR_NUM_IN);
    int_in[INTEGRATOR_P] = U[k];
    int_in[INTEGRATOR_X0] = X[k];

    // Create an evaluation node
    vector<MX> I_out = integrator.call(int_in);

    // Save continuity constraints
    g[k] = I_out[INTEGRATOR_XF] - X[k+1];
  }
  
  // NLP objective function
  MX ff = mterm.call(X.back()).front();
  MXFunction F(V,ff);

  // NLP constraint function
  MX gg = vertcat(g);
  MXFunction G(V,gg);
  G.setOption("number_of_fwd_dir",7);
  G.setOption("ad_mode","forward");
  G.setOption("numeric_jacobian",false);
  G.setOption("name","NLP constraint function");
  G.init();
  
  FX h_jac;
  if(exact_hessian){
    // Lagrange multiplier
    MX lag("lag",gg.size1());

    // Objective function scaling factor
    MX sigma("sigma");
    
    // Lagrangian function
    vector<MX> lfcn_in(3);
    lfcn_in[0] = V;
    lfcn_in[1] = lag;
    lfcn_in[2] = sigma;
    MXFunction lfcn(lfcn_in,sigma*ff + inner_prod(lag,gg));
    lfcn.init();
    
    // Gradient of the lagrangian
    vector<MX> lgrad = lfcn.grad();
    MXFunction lgfcn(lfcn_in,trans(lgrad[0]));
    lgfcn.init();

    // Hessian of the lagrangian
    h_jac = lgfcn.jacobian();
    h_jac.init();
  }
  
  // Create an NLP solver instance
  IpoptSolver nlp_solver(F,G,h_jac);
  nlp_solver.setOption("tol",1e-5);
  if(!exact_hessian)
    nlp_solver.setOption("hessian_approximation", "limited-memory");
  nlp_solver.setOption("max_iter",100);
  nlp_solver.setOption("linear_solver","ma57");
  // nlp_solver.setOption("verbose",true);
  // nlp_solver.setOption("derivative_test","first-order");
  nlp_solver.init();
    
  // Initial guess and bounds on variables
  Matrix<double>& V_init = nlp_solver.input(NLP_X_INIT);
  Matrix<double>& V_min  = nlp_solver.input(NLP_LBX);
  Matrix<double>& V_max  = nlp_solver.input(NLP_UBX);
  
  // All nonlinear constraints are equality constraints
  nlp_solver.input(NLP_LBG).setAll(0);
  nlp_solver.input(NLP_UBG).setAll(0);
  
  // Running index
  el=0;
  
  for(int k=0; k<=ns; ++k){
    // Pass bounds and guess for state
    for(int i=0; i<nx; ++i){
      V_init.at(el) = x_init[i];
      V_min.at(el) = k == 0 ? x0_min[i] : k==ns ? xf_min[i] : x_min[i];
      V_max.at(el) = k == 0 ? x0_max[i] : k==ns ? xf_max[i] : x_max[i];
      el++;
    }
    
    // Variables below do not appear at the end point
    if(k==ns) break;

    // Pass bounds and guess for control
    for(int i=0; i<nu; ++i){
      V_init.at(el) = u_init[i];
      V_min.at(el) = u_min[i];
      V_max.at(el) = u_max[i];
      el++;
    }
  }
  
  // Make sure all values set
  casadi_assert(el==NV);
  
  // Solve the problem
  nlp_solver.solve();

  // Optimal solution of the NLP
  const Matrix<double>& V_opt = nlp_solver.output(NLP_X_OPT);
  
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
