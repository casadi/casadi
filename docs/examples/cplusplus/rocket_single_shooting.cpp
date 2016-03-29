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
extern "C" void casadi_load_nlpsol_ipopt();
extern "C" void casadi_load_nlpsol_scpgen();

bool sundials_integrator = true;
bool explicit_integrator = false;
bool lifted_newton = false;

int main(){
  // Load integrators manually
  casadi_load_integrator_cvodes();
  casadi_load_integrator_idas();
  casadi_load_integrator_rk();
  casadi_load_nlpsol_ipopt();
  casadi_load_nlpsol_scpgen();

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
  SX x = SX::vertcat({s, v, m});

  // Control
  SX u = SX::sym("u");

  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Differential equation
  SX rhs = SX::vertcat({v, (u-alpha*v*v)/m, -beta*u*u});

  // Initial conditions
  vector<double> x0 = {0, 0, 1};

  // DAE
  SXDict dae = {{"x", x}, {"p", u}, {"t", t}, {"ode", rhs}};

  // Integrator options
  string plugin;
  Dict opts;
  if(sundials_integrator){
    if(explicit_integrator){
      // Explicit integrator (CVODES)
      plugin = "cvodes";
      // opts["exact_jacobian"] = true;
      // opts["linear_multistep_method"] = "bdf"; // adams or bdf
      // opts["nonlinear_solver_iteration"] = "newton"; // newton or functional
    } else {
      // Implicit integrator (IDAS)
      plugin = "idas";
      opts["calc_ic"] = false;
    }
    opts["fsens_err_con"] = true;
    opts["quad_err_con"] = true;
    opts["abstol"] = 1e-6;
    opts["reltol"] = 1e-6;
    opts["stop_at_end"] = false;
    //  opts["fsens_all_at_once"] = false;
    opts["steps_per_checkpoint"] = 100; // BUG: Too low number causes segfaults
  } else {
    // An explicit Euler integrator
    plugin = "rk";
    opts["expand_f"] = true;
    opts["interpolation_order"] = 1;
    opts["number_of_finite_elements"] = 1000;
  }
  opts["t0"] = t0;
  opts["tf"] = tf;

  // Create integrator
  Function F = integrator("integrator", plugin, dae, opts);

  // control for all segments
  MX U = MX::sym("U",nu);

  // Integrate over all intervals
  MX X=X0;
  for(int k=0; k<nu; ++k){
    // Integrate
    X = F(MXDict{{"x0", X}, {"p", U(k)}}).at("xf");

    // Lift X
    if(lifted_newton){
      X = lift(X, X);
    }
  }

  // Objective function
  MX J = dot(U,U);

  // Terminal constraints
  MX G = vertcat(X(0),X(1));

  // Create the NLP
  MXDict nlp = {{"x", U}, {"f", J}, {"g", G}};

  // NLP solver options
  Dict solver_opts;
  string solver_name;
  if(lifted_newton){
    solver_name = "scpgen";
    solver_opts["verbose"] = true;
    solver_opts["regularize"] = false;
    solver_opts["max_iter_ls"] = 1;
    solver_opts["max_iter"] = 100;
    solver_opts["qpsol"] = "nlp"; // Use IPOPT as QP solver
    Dict ipopt_options;
    ipopt_options["tol"] = 1e-12;
    ipopt_options["print_level"] = 0;
    ipopt_options["print_time"] = false;
    solver_opts["qpsol_options"] =
      Dict{{"nlpsol_options", Dict{{"nlpsol", "ipopt"}}}};
  } else {
    solver_name = "ipopt";
    solver_opts["ipopt.tol"] = 1e-10;
    solver_opts["ipopt.hessian_approximation"] = "limited-memory";
  }

  // Create NLP solver
  Function solver = nlpsol("nlpsol", solver_name, nlp, solver_opts);

  // Bounds on u and initial condition
  vector<double> umin(nu, -10), umax(nu, 10), u0(nu, 0.4);

  // Bounds on g
  vector<double> gmin = {10, 0}, gmax = {10, 0};

  // Solve NLP
  vector<double> Usol;
  solver({{"lbx", umin}, {"ubx", umax}, {"x0", u0}, {"lbg", gmin}, {"ubg", gmax}},
         {{"x", &Usol}});
  cout << "optimal solution: " << Usol << endl;

  return 0;
}
