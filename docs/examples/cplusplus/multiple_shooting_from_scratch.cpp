/*
 *    MIT No Attribution
 *
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
#include <casadi/casadi.hpp>

using namespace casadi;

int main(){
  // Declare variables
  SX u = SX::sym("u"); // control
  SX r = SX::sym("r"), s = SX::sym("s"); // states
  SX x = vertcat(r,s);

  // Number of differential states
  int nx = x.size1();

  // Number of controls
  int nu = u.size1();

  // Bounds and initial guess for the control
  std::vector<double> u_min =  { -0.75 };
  std::vector<double> u_max  = {  1.0  };
  std::vector<double> u_init = {  0.0  };

  // Bounds and initial guess for the state
  std::vector<double> x0_min = {   0,    1 };
  std::vector<double> x0_max = {   0,    1 };
  std::vector<double> x_min  = {-inf, -inf };
  std::vector<double> x_max  = { inf,  inf };
  std::vector<double> xf_min = {   0,    0 };
  std::vector<double> xf_max = {   0,    0 };
  std::vector<double> x_init = {   0,    0 };

  // Final time
  double tf = 20.0;

  // Number of shooting nodes
  int ns = 50;

  // ODE right hand side and quadrature
  SX ode = vertcat((1 - s*s)*r - s + u, r);
  SX quad = r*r + s*s + u*u;
  SXDict dae = {{"x", x}, {"p", u}, {"ode", ode}, {"quad", quad}};

  // Create an integrator (CVodes)
  Function F = integrator("integrator", "cvodes", dae, 0, tf/ns);

  // Total number of NLP variables
  int NV = nx*(ns+1) + nu*ns;

  // Declare variable vector for the NLP
  MX V = MX::sym("V",NV);

  // NLP variable bounds and initial guess
  std::vector<double> v_min,v_max,v_init;

  // Offset in V
  int offset=0;

  // State at each shooting node and control for each shooting interval
  std::vector<MX> X, U;
  for(int k=0; k<ns; ++k){
    // Local state
    X.push_back( V.nz(Slice(offset,offset+nx)));
    if(k==0){
      v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
      v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
    } else {
      v_min.insert(v_min.end(), x_min.begin(), x_min.end());
      v_max.insert(v_max.end(), x_max.begin(), x_max.end());
    }
    v_init.insert(v_init.end(), x_init.begin(), x_init.end());
    offset += nx;

    // Local control
    U.push_back( V.nz(Slice(offset,offset+nu)));
    v_min.insert(v_min.end(), u_min.begin(), u_min.end());
    v_max.insert(v_max.end(), u_max.begin(), u_max.end());
    v_init.insert(v_init.end(), u_init.begin(), u_init.end());
    offset += nu;
  }

  // State at end
  X.push_back(V.nz(Slice(offset,offset+nx)));
  v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
  v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
  v_init.insert(v_init.end(), x_init.begin(), x_init.end());
  offset += nx;

  // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
  casadi_assert(offset==NV, "");

  // Objective function
  MX J = 0;

  //Constraint function and bounds
  std::vector<MX> g;

  // Loop over shooting nodes
  for(int k=0; k<ns; ++k){
    // Create an evaluation node
    MXDict I_out = F(MXDict{{"x0", X[k]}, {"p", U[k]}});

    // Save continuity constraints
    g.push_back( I_out.at("xf") - X[k+1] );

    // Add objective function contribution
    J += I_out.at("qf");
  }

  // NLP
  MXDict nlp = {{"x", V}, {"f", J}, {"g", vertcat(g)}};

  // Set options
  Dict opts;
  opts["ipopt.tol"] = 1e-5;
  opts["ipopt.max_iter"] = 100;

  // Create an NLP solver and buffers
  Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);
  std::map<std::string, DM> arg, res;

  // Bounds and initial guess
  arg["lbx"] = v_min;
  arg["ubx"] = v_max;
  arg["lbg"] = 0;
  arg["ubg"] = 0;
  arg["x0"] = v_init;

  // Solve the problem
  res = solver(arg);

  // Optimal solution of the NLP
  std::vector<double> V_opt(res.at("x"));

  // Get the optimal state trajectory
  std::vector<double> r_opt(ns+1), s_opt(ns+1);
  for(int i=0; i<=ns; ++i){
    r_opt[i] = V_opt.at(i*(nx+1));
    s_opt[i] = V_opt.at(1+i*(nx+1));
  }
  std::cout << "r_opt = " << std::endl << r_opt << std::endl;
  std::cout << "s_opt = " << std::endl << s_opt << std::endl;

  // Get the optimal control
  std::vector<double> u_opt(ns);
  for(int i=0; i<ns; ++i){
    u_opt[i] = V_opt.at(nx + i*(nx+1));
  }
  std::cout << "u_opt = " << std::endl << u_opt << std::endl;


  return 0;
}
