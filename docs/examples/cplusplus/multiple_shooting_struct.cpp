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


/**
  The problem of multiple_shooting_from_scratch.cpp, with the NLP variables, bounds and
  constraints laid out by a Struct instead of manual offset bookkeeping.
*/

#include <casadi/casadi.hpp>

using namespace casadi;

int main() {
  // ODE with states r, s and control u
  SX r = SX::sym("r"), s = SX::sym("s"), u = SX::sym("u");
  SXDict dae = {{"x", vertcat(r, s)}, {"p", u},
                {"ode", vertcat((1 - s*s)*r - s + u, r)}, {"quad", r*r + s*s + u*u}};
  double tf = 20.0;
  casadi_int ns = 50;
  Function F = integrator("integrator", "cvodes", dae, 0, tf/ns);

  // NLP variables: states at ns+1 nodes and controls on ns intervals, interleaved
  Struct V;
  V.add("X", 2, 1, {ns+1});
  V.add("U", 1, 1, {ns});
  V.interleave({"X", "U"});

  // Symbolic variables, bounds and initial guess
  StructMX v = StructMX::sym(V);
  StructDM lbx(V, -inf), ubx(V, inf), x0(V, 0);
  lbx("U") = -0.75;
  ubx("U") = 1.0;
  lbx("X", 0) = ubx("X", 0) = DM({0, 1});
  lbx("X", -1) = ubx("X", -1) = DM({0, 0});

  // Continuity constraints and objective
  Struct G;
  G.add("cont", 2, 1, {ns});
  StructMX g(G, 0);
  MX J = 0;
  for (casadi_int k=0; k<ns; ++k) {
    MXDict I = F(MXDict{{"x0", v("X", k)}, {"p", v("U", k)}});
    g("cont", k) = I.at("xf") - v("X", k+1);
    J += I.at("qf");
  }

  // Solve
  Function solver = nlpsol("nlpsol", "ipopt", {{"x", v.cat()}, {"f", J}, {"g", g.cat()}},
    Dict{{"ipopt.tol", 1e-5}, {"ipopt.max_iter", 100}});
  DMDict res = solver(DMDict{{"lbx", lbx.cat()}, {"ubx", ubx.cat()}, {"x0", x0.cat()},
    {"lbg", 0}, {"ubg", 0}});

  // Trajectories, one column per node
  StructDM sol(V, res.at("x"));
  std::cout << "r_opt = " << sol("X", ":", 0) << std::endl;
  std::cout << "s_opt = " << sol("X", ":", 1) << std::endl;
  std::cout << "u_opt = " << sol("U") << std::endl;
  return 0;
}
