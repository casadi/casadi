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
   Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
   @author: Joel Andersson, KU Leuven 2013
*/

#include <casadi/casadi.hpp>
#include <iomanip>

using namespace casadi;

int main(int argc, char *argv[]) {
  // End time
  double tf = 10.0;

  // Dimensions
  int nx = 3;
  int np = 1;

  // Declare variables
  SX x  = SX::sym("x",nx);  // state
  SX p  = SX::sym("u",np);  // control

  // ODE right hand side function
  SX ode = vertcat((1 - x(1)*x(1))*x(0) - x(1) + p,
                   x(0),
                   x(0)*x(0) + x(1)*x(1) + p*p);
  SXDict dae = {{"x", x}, {"p", p}, {"ode", ode}};

  // Number of finite elements
  int n = 100;

  // Size of the finite elements
  double h = tf/n;

  // Degree of interpolating polynomial
  int d = 4;

  // Choose collocation points
  std::vector<double> tau_root = collocation_points(d, "legendre");
  tau_root.insert(tau_root.begin(), 0);

  // Nonlinear solver to use
  std::string solver = "newton";
  if (argc>1) solver = argv[1]; // chose a different solver from command line

  // Coefficients of the collocation equation
  std::vector<std::vector<double> > C(d+1, std::vector<double>(d+1,0));

  // Coefficients of the continuity equation
  std::vector<double> D(d+1,0);

  // For all collocation points
  for(int j=0; j<d+1; ++j){

    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    Polynomial p = 1;
    for(int r=0; r<d+1; ++r){
      if(r!=j){
        p *= Polynomial(-tau_root[r],1)/(tau_root[j]-tau_root[r]);
      }
    }

    // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D[j] = p(1.0);

    // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    Polynomial dp = p.derivative();
    for(int r=0; r<d+1; ++r){
      C[j][r] = dp(tau_root[r]);
    }
  }

  // Total number of variables for one finite element
  MX X0 = MX::sym("X0",nx);
  MX P  = MX::sym("P",np);
  MX V = MX::sym("V",d*nx);

  // Get the state at each collocation point
  std::vector<MX> X(1,X0);
  for(int r=0; r<d; ++r){
    X.push_back(V.nz(Slice(r*nx,(r+1)*nx)));
  }

  // Get the collocation equations (that define V)
  Function f("f", {dae["x"], dae["p"]}, {dae["ode"]});
  std::vector<MX> V_eq;
  for(int j=1; j<d+1; ++j){
    // Expression for the state derivative at the collocation point
    MX xp_j = 0;
    for(int r=0; r<d+1; ++r){
      xp_j += C[r][j]*X[r];
    }

    // Append collocation equations
    std::vector<MX> v = {X[j], P};
    v = f(v);
    V_eq.push_back(h*v[0] - xp_j);
  }

  // Root-finding function, implicitly defines V as a function of X0 and P
  Function vfcn("vfcn", {V, X0, P}, {vertcat(V_eq)});

  // Convert to sxfunction to decrease overhead
  Function vfcn_sx = vfcn.expand("vfcn");

  // Create a implicit function instance to solve the system of equations
  Dict opts;
  if (solver=="ipopt") {
    // Use an NLP solver
    opts["nlpsol"] = "ipopt";
    opts["nlpsol_options"] = Dict{{"print_time", false}, {"ipopt.print_level", 0}};
    solver = "nlpsol";
  } else if (solver=="kinsol") {
    opts["linear_solver_type"] = "user_defined";
  }
  Function ifcn = rootfinder("ifcn", solver, vfcn_sx, opts);

  // Get an expression for the state at the end of the finite element
  std::vector<MX> ifcn_arg = {MX(), X0, P};
  V = ifcn(ifcn_arg).front();
  X.resize(1);
  for(int r=0; r<d; ++r){
    X.push_back(V.nz(Slice(r*nx, (r+1)*nx)));
  }
  MX XF = 0;
  for(int r=0; r<d+1; ++r){
    XF += D[r]*X[r];
  }

  // Get the discrete time dynamics
  Function F("F", {X0, P}, {XF});

  // Do this iteratively for all finite elements
  MX Xk = X0;
  for(int i=0; i<n; ++i){
    Xk = F(std::vector<MX>{Xk, P}).at(0);
  }

  // Fixed-step integrator
  Function irk_integrator("irk_integrator", MXDict{{"x0", X0}, {"p", P}, {"xf", Xk}},
                          integrator_in(), integrator_out());

  // Create a conventional integrator for reference
  Function ref_integrator = integrator("ref_integrator",
                                       "cvodes", dae, {{"tf", tf}});

  // Test values
  std::vector<double> x0_val = {0, 1, 0};
  double p_val = 0.2;

  // Make sure that both integrators give consistent results
  for(int integ=0; integ<2; ++integ){
    Function F = integ==0 ? irk_integrator : ref_integrator;
    std::cout << "-------" << std::endl;
    std::cout << "Testing " << F.name() << std::endl;
    std::cout << "-------" << std::endl;

    // Generate a new function that calculates forward and reverse directional derivatives
    Function dF = F.factory("dF", {"x0", "p", "fwd:x0", "fwd:p", "adj:xf"},
                                  {"xf", "fwd:xf", "adj:x0", "adj:p"});

    // Arguments for evaluation
    std::map<std::string, DM> arg, res;
    arg["x0"] = x0_val;
    arg["p"] = p_val;

    // Forward sensitivity analysis, first direction: seed p and x0[0]
    arg["fwd_x0"] = std::vector<double>{1, 0, 0};
    arg["fwd_p"] = 1;

    // Adjoint sensitivity analysis, seed xf[2]
    arg["adj_xf"] = std::vector<double>{0, 0, 1};

    // Integrate
    res = dF(arg);

    // Get the nondifferentiated results
    std::cout << std::setw(15) << "xf = " << res.at("xf") << std::endl;

    // Get the forward sensitivities
    std::cout << std::setw(15) << "d(xf)/d(p)+d(xf)/d(x0[0]) = " <<  res.at("fwd_xf") << std::endl;

    // Get the adjoint sensitivities
    std::cout << std::setw(15) << "d(xf[2])/d(x0) = " << res.at("adj_x0") << std::endl;
    std::cout << std::setw(15) << "d(xf[2])/d(p) = " << res.at("adj_p") << std::endl;
  }
  return 0;
}
