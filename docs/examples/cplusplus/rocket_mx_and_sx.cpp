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

// CONSTRUCT THE INTEGRATOR
Function create_integrator(int nj, int nu){
  SX u = SX::sym("u"); // control for one segment

  // Initial position
  SX s0 = SX::sym("s0"); // initial position
  SX v0 = SX::sym("v0"); // initial speed
  SX m0 = SX::sym("m0"); // initial mass

  SX dt = 10.0/(nj*nu); // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Integrate over the interval with Euler forward
  SX s = s0, v = v0, m = m0;

  SX dm = -dt*beta*u*u;
  for(int j=0; j<nj; ++j){
    s += dt*v;
    v += dt / m * (u - alpha * v*v);
    m += dm;
  }

  // State vector
  SX x = vertcat(s, v, m);
  SX x0 = vertcat(s0, v0, m0);

  // Integrator
  return Function("integrator", {u, x0}, {x});
}


int main(){
  // Dimensions
  int nj = 1000; // Number of integration steps per control segment
  int nu = 1000; // Number of control segments

  // Create an integrator
  Function integrator = create_integrator(nj,nu);

  // PART 2: CONSTRUCT THE NLP
  MX U = MX::sym("U",nu); // control for all segments

  // Initial position
  vector<double> X0(3);
  X0[0] = 0; // initial position
  X0[1] = 0; // initial speed
  X0[2] = 1; // initial mass

  // Integrate over all intervals
  MX X=X0;
  for(int k=0; k<nu; ++k){
    // Assemble the input
    vector<MX> input(2);
    input[0] = U(k);
    input[1] = X;

    // Integrate
    X = integrator(input).at(0);
  }

  // Objective function
  MX F = dot(U,U);

  // Terminal constraints
  MX G = vertcat(X(0),X(1));

  // Create the NLP
  MXDict nlp = {{"x", U}, {"f", F}, {"g", G}};

  // Allocate an NLP solver and buffers
  Dict opts = {{"ipopt.tol", 1e-10},
               {"ipopt.hessian_approximation", "limited-memory"}};
  Function solver = nlpsol("solver", "ipopt", nlp, opts);
  std::map<std::string, DM> arg, res;

  // Bounds on u and initial condition
  arg["lbx"] = -10;
  arg["ubx"] = 10;
  arg["x0"] = 0.4;

  // Bounds on g
  vector<double> Gmin(2), Gmax(2);
  Gmin[0] = Gmax[0] = 10;
  Gmin[1] = Gmax[1] =  0;
  arg["lbg"] = Gmin;
  arg["ubg"] = Gmax;

  // Solve the problem
  res = solver(arg);

  // Get the solution
  cout << "optimal solution: " << res.at("x") << endl;

  return 0;
}
