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


#include <casadi/casadi.hpp>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace casadi;

struct Test {
  SXDict dae;
  double tf;
  vector<double> x0;
  double u0;
  bool is_ode;
  string name;
};

/** \brief Generate a simple ODE */
Test simpleODE(){
  // Time
  SX t = SX::sym("t");

  // Parameter
  SX u = SX::sym("u");

  // Differential states
  SX s = SX::sym("s"), v = SX::sym("v"), m = SX::sym("m");
  SX x = vertcat(s, v, m);

  // Constants
  double alpha = 0.05; // friction
  double beta = 0.1;   // fuel consumption rate

  // Differential equation
  SX ode = vertcat(v, (u-alpha*v*v)/m, -beta*u*u);

  // Quadrature
  SX quad = pow(v,3) + pow((3-sin(t))-u,2);

  // Return test problem
  Test r;
  r.dae = decltype(r.dae){{"t", t}, {"x", x}, {"p", u}, {"ode", ode}, {"quad", quad}}; // decltype because MSVC bug
  r.tf = 0.5;
  r.x0 = {0, 0, 1};
  r.u0 = 0.4;
  r.is_ode = true;
  r.name = "simple_ode";
  return r;
}

/** \brief Generate a simple DAE */
Test simpleDAE(){
  // Parameter
  SX u = SX::sym("u");

  // Differential state
  SX x = SX::sym("x");

  // Algebraic variable
  SX z = SX::sym("z");

  // Differential equation
  SX ode = -x + 0.5*x*x + u + 0.5*z;

  // Algebraic constraint
  SX alg = z + exp(z) - 1. + x;

  // Quadrature
  SX quad = x*x + 3.0*u*u;

  // Return DAE
  Test r;
  r.dae = decltype(r.dae){{"x", x}, {"z", z}, {"p", u}, {"ode", ode}, {"alg", alg}, {"quad", quad}};
  r.tf = 5;
  r.x0 = {1};
  r.u0 = 0.4;
  r.is_ode = false;
  r.name = "simple_dae";
  return r;
}

struct Solver {
  string plugin;
  bool ode_only;
  Dict opts;
};

int main(){

  // Test problems
  vector<Test> tests = {simpleODE(), simpleDAE()};

  // ODE/DAE integrators
  vector<Solver> solvers;
  solvers.push_back({"cvodes", true, Dict()});
  solvers.push_back({"idas", false, Dict()});
  solvers.push_back({"rk", true, Dict()});
  Dict kinsol_options = {{"linear_solver_type", "dense"}};

  Dict coll_opts = {{"rootfinder", "kinsol"},
                    {"collocation_scheme", "legendre"},
                    {"rootfinder_options", kinsol_options}};
  solvers.push_back({"collocation", false, coll_opts});

  // Loop over all problems
  for (auto&& test : tests) {
    // Loop over all solvers
    for (auto&& solver : solvers) {
      // Skip if problem cannot be handled
      if (solver.ode_only && !test.is_ode) continue;
      // Printout
      cout << "Solving \"" << test.name << "\" using \"" << solver.plugin << "\"" << endl;

      // Get integrator
      Dict opts = solver.opts;
      opts["tf"] = test.tf;
      Function I = integrator("I", solver.plugin, test.dae, opts);

      // Buffers for evaluation
      std::map<std::string, DM> arg, res;

      // Integrate to get results
      arg = decltype(arg){{"x0", test.x0},
                          {"p", test.u0}};
      res = I(arg);
      vector<double> xf(res.at("xf"));
      vector<double> qf(res.at("qf"));
      cout << setw(50) << "Unperturbed solution: " << "xf  = " << xf <<  ", qf  = " << qf << endl;

      // Perturb solution to get a finite difference approximation
      double h = 0.001;
      arg["p"] = test.u0+h;
      res = I(arg);
      vector<double> fd_xf((res.at("xf")-xf)/h);
      vector<double> fd_qf((res.at("qf")-qf)/h);
      cout << setw(50) << "Finite difference approximation: " << "d(xf)/d(p) = " << fd_xf << ", d(qf)/d(p) = " << fd_qf << endl;

      // Calculate once, forward
      Function I_fwd = I.factory("I_fwd", {"x0", "p", "fwd:x0", "fwd:p"},
                                          {"fwd:xf", "fwd:qf"});
      arg = decltype(arg){{"x0", test.x0},
                          {"p", test.u0},
                          {"fwd_x0", 0},
                          {"fwd_p", 1}};
      res = I_fwd(arg);
      vector<double> fwd_xf(res.at("fwd_xf"));
      vector<double> fwd_qf(res.at("fwd_qf"));
      cout << setw(50) << "Forward sensitivities: " << "d(xf)/d(p) = " << fwd_xf << ", d(qf)/d(p) = " << fwd_qf << endl;

      // Calculate once, adjoint
      Function I_adj = I.factory("I_adj", {"x0", "p", "adj:xf", "adj:qf"},
                                          {"adj:p", "adj:x0"});
      arg = decltype(arg){{"x0", test.x0}, {"p", test.u0}, {"adj_xf", 0}, {"adj_qf", 1}};
      res = I_adj(arg);
      vector<double> adj_x0(res.at("adj_x0"));
      vector<double> adj_p(res.at("adj_p"));
      cout << setw(50) << "Adjoint sensitivities: " << "d(qf)/d(x0) = " << adj_x0 << ", d(qf)/d(p) = " << adj_p << endl;

      // Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
      arg["p"] = test.u0+h;
      res = I_adj(arg);
      vector<double> fd_adj_x0((res.at("adj_x0")-adj_x0)/h);
      vector<double> fd_adj_p((res.at("adj_p")-adj_p)/h);
      cout << setw(50) << "FD of adjoint sensitivities: " << "d2(qf)/d(x0)d(p) = " << fd_adj_x0 << ", d2(qf)/d(p)d(p) = " << fd_adj_p << endl;

      // Forward over adjoint to get the second order sensitivities
      Function I_foa = I_adj.factory("I_foa",
                                     {"x0", "p", "fwd:p", "adj_xf", "adj_qf"},
                                     {"fwd:adj_x0", "fwd:adj_p"});
      arg = decltype(arg){{"x0", test.x0},
                          {"p", test.u0},
                          {"fwd_p", 1},
                          {"adj_xf", 0},
                          {"adj_qf", 1}};
      res = I_foa(arg);
      vector<double> fwd_adj_x0(res.at("fwd_adj_x0"));
      vector<double> fwd_adj_p(res.at("fwd_adj_p"));
      cout << setw(50) << "Forward over adjoint sensitivities: " << "d2(qf)/d(x0)d(p) = " << fwd_adj_x0 << ", d2(qf)/d(p)d(p) = " << fwd_adj_p << endl;

      // Adjoint over adjoint to get the second order sensitivities
      Function I_aoa = I_adj.factory("I_aoa", {"x0", "p", "adj_xf", "adj_qf", "adj:adj_p"},
                                              {"adj:x0", "adj:p"});
      arg = decltype(arg){{"x0", test.x0},
                          {"p", test.u0},
                          {"adj_xf", 0},
                          {"adj_qf", 1},
                          {"adj_adj_p", 1}};
      res = I_aoa(arg);
      vector<double> adj_adj_x0(res.at("adj_x0"));
      vector<double> adj_adj_p(res.at("adj_p"));
      cout << setw(50) << "Adjoint over adjoint sensitivities: " << "d2(qf)/d(x0)d(p) = " << adj_adj_x0 << ", d2(qf)/d(p)d(p) = " << adj_adj_p << endl;
    }
  }
  return 0;
}
