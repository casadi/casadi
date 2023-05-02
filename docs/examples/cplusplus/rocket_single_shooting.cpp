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

#include <iostream>
#include <casadi/casadi.hpp>

using namespace casadi;

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
  std::vector<double> X0(3);
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
  std::vector<double> x0 = {0, 0, 1};

  // DAE
  SXDict dae = {{"x", x}, {"p", u}, {"t", t}, {"ode", rhs}};

  // Integrator options
  std::string plugin;
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

  // Create integrator
  Function F = integrator("integrator", plugin, dae, t0, tf, opts);

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
  std::string solver_name;
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
  std::vector<double> umin(nu, -10), umax(nu, 10), u0(nu, 0.4);

  // Bounds on g
  std::vector<double> gmin = {10, 0}, gmax = {10, 0};

  // Solve NLP
  std::vector<double> Usol;
  solver({{"lbx", umin}, {"ubx", umax}, {"x0", u0}, {"lbg", gmin}, {"ubg", gmax}},
         {{"x", &Usol}});
  std::cout << "optimal solution: " << Usol << std::endl;

  return 0;
}
