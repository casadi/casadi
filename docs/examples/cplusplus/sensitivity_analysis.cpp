/*
 *    MIT No Attribution
 *
 *    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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

#include <casadi/casadi.hpp>

#include <iostream>
#include <iomanip>

using namespace casadi;

struct Test {
  SXDict dae;
  double tf;
  std::vector<double> x0;
  double u0;
  bool is_ode;
  std::string name;
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
  std::string plugin;
  bool ode_only;
  Dict opts;
};

int main(){

  // Test problems
  std::vector<Test> tests = {simpleODE(), simpleDAE()};

  // ODE/DAE integrators
  std::vector<Solver> solvers;
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
      // Test with or without multiple output times
      for (int ntout : {1, 3}) {
        // Skip if problem cannot be handled
        if (solver.ode_only && !test.is_ode) continue;
        // Printout
        std::cout << "Solving \"" << test.name << "\" using \"" << solver.plugin << "\", " << ntout << " output times" << std::endl;
        // Output time grid
        std::vector<double> tout;
        if (ntout == 1) {
          tout = {test.tf};
        } else {
          tout = {0, 0.5 * test.tf, test.tf};
        }
        // Create integrator instance
        Function I = integrator("I", solver.plugin, test.dae, 0, tout, solver.opts);

        // Buffers for evaluation
        std::map<std::string, DM> arg, res;

        // Integrate to get results
        arg = decltype(arg){{"x0", test.x0},
                            {"p", test.u0}};
        res = I(arg);
        DM xf(res.at("xf"));
        DM qf(res.at("qf"));
        std::cout << std::setw(50) << "Unperturbed solution: " << "xf  = " << xf.nonzeros() <<  ", qf  = " << qf.nonzeros() << std::endl;

        // Perturb solution to get a finite difference approximation
        double h = 0.001;
        arg["p"] = test.u0+h;
        res = I(arg);
        DM fd_xf = (res.at("xf")-xf)/h;
        DM fd_qf = (res.at("qf")-qf)/h;
        std::cout << std::setw(50) << "Finite difference approximation: " << "d(xf)/d(p) = " << fd_xf.nonzeros() << ", d(qf)/d(p) = " << fd_qf.nonzeros() << std::endl;

        // Calculate once, forward
        Function I_fwd = I.factory("I_fwd", {"x0", "p", "fwd:x0", "fwd:p"},
                                            {"fwd:xf", "fwd:qf"});
        arg = decltype(arg){{"x0", test.x0},
                            {"p", test.u0},
                            {"fwd_x0", 0},
                            {"fwd_p", 1}};
        res = I_fwd(arg);
        DM fwd_xf = res.at("fwd_xf");
        DM fwd_qf = res.at("fwd_qf");
        std::cout << std::setw(50) << "Forward sensitivities: " << "d(xf)/d(p) = " << fwd_xf.nonzeros() << ", d(qf)/d(p) = " << fwd_qf.nonzeros() << std::endl;

        // Calculate once, adjoint
        std::vector<double> adj_qf(ntout, 0);
        adj_qf.back() = 1;
        Function I_adj = I.factory("I_adj", {"x0", "p", "adj:xf", "adj:qf"},
                                            {"adj:p", "adj:x0"});
        arg = decltype(arg){{"x0", test.x0}, {"p", test.u0}, {"adj_xf", 0}, {"adj_qf", adj_qf}};
        res = I_adj(arg);
        DM adj_x0 = res.at("adj_x0");
        DM adj_p = res.at("adj_p");
        std::cout << std::setw(50) << "Adjoint sensitivities: " << "d(qf)/d(x0) = " << adj_x0.nonzeros() << ", d(qf)/d(p) = " << adj_p.nonzeros() << std::endl;

        // Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
        arg["p"] = test.u0+h;
        res = I_adj(arg);
        DM fd_adj_x0 = (res.at("adj_x0")-adj_x0)/h;
        DM fd_adj_p = (res.at("adj_p")-adj_p)/h;
        std::cout << std::setw(50) << "FD of adjoint sensitivities: " << "d2(qf)/d(x0)d(p) = " << fd_adj_x0.nonzeros() << ", d2(qf)/d(p)d(p) = " << fd_adj_p.nonzeros() << std::endl;

        // Forward over adjoint to get the second order sensitivities
        Function I_foa = I_adj.factory("I_foa",
                                      {"x0", "p", "fwd:p", "adj_xf", "adj_qf"},
                                      {"fwd:adj_x0", "fwd:adj_p"});
        arg = decltype(arg){{"x0", test.x0},
                            {"p", test.u0},
                            {"fwd_p", 1},
                            {"adj_xf", 0},
                            {"adj_qf", adj_qf}};
        res = I_foa(arg);
        DM fwd_adj_x0 = res.at("fwd_adj_x0");
        DM fwd_adj_p = res.at("fwd_adj_p");
        std::cout << std::setw(50) << "Forward over adjoint sensitivities: " << "d2(qf)/d(x0)d(p) = " << fwd_adj_x0.nonzeros() << ", d2(qf)/d(p)d(p) = " << fwd_adj_p.nonzeros() << std::endl;

        // Adjoint over adjoint to get the second order sensitivities
        Function I_aoa = I_adj.factory("I_aoa", {"x0", "p", "adj_xf", "adj_qf", "adj:adj_p"},
                                                {"adj:x0", "adj:p"});
        arg = decltype(arg){{"x0", test.x0},
                            {"p", test.u0},
                            {"adj_xf", 0},
                            {"adj_qf", adj_qf},
                            {"adj_adj_p", 1}};
        res = I_aoa(arg);
        DM adj_adj_x0 = res.at("adj_x0");
        DM adj_adj_p = res.at("adj_p");
        std::cout << std::setw(50) << "Adjoint over adjoint sensitivities: " << "d2(qf)/d(x0)d(p) = " << adj_adj_x0.nonzeros() << ", d2(qf)/d(p)d(p) = " << adj_adj_p.nonzeros() << std::endl;
      }
    }
  }
  return 0;
}
