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

#include <iostream>
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/fx/external_function.hpp>

using namespace CasADi;
using namespace std;

// CONSTRUCT THE INTEGRATOR
FX create_integrator(int nj, int nu){
  SX u("u"); // control for one segment

  // Initial position
  SX s0("s0"); // initial position
  SX v0("v0"); // initial speed
  SX m0("m0"); // initial mass

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
  SXMatrix x, x0;
  x.append(s);
  x.append(v);
  x.append(m);
  x0.append(s0);
  x0.append(v0);
  x0.append(m0);

  // Integrator
  vector<SXMatrix> input(2);
  input[0] = u;
  input[1] = x0;
  SXMatrix output = x;
  SXFunction integrator(input,output);
  integrator.init();

  return integrator;
}


int main(){
  {
  // Dimensions
  int nj = 1000; // Number of integration steps per control segment
  int nu = 1000; // Number of control segments

  // Create an integrator
  FX integrator = create_integrator(nj,nu);

  // PART 2: CONSTRUCT THE NLP
  MX U("U",nu); // control for all segments
 
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
    input[0] = U[k];
    input[1] = X;

    // Integrate
    X = integrator.call(input).at(0);
  }

  // Objective function
  MX F = inner_prod(U,U);

  // Terminal constraints
  MX G = vertcat(X[0],X[1]);
  
  // Create the NLP
  MXFunction nlp(nlpIn("x",U),nlpOut("f",F,"g",G));

  // Allocate an NLP solver
  IpoptSolver solver(nlp);
  
  // Set options
  solver.setOption("tol",1e-10);
  solver.setOption("hessian_approximation","limited-memory");

  // initialize the solver
  solver.init();

  // Bounds on u and initial condition
  vector<double> Umin(nu), Umax(nu), Usol(nu);
  for(int i=0; i<nu; ++i){
    Umin[i] = -10;
    Umax[i] =  10;
    Usol[i] = 0.4;
  }
  solver.setInput(Umin,"lbx");
  solver.setInput(Umax,"ubx");
  solver.setInput(Usol,"x0");

  // Bounds on g
  vector<double> Gmin(2), Gmax(2);
  Gmin[0] = Gmax[0] = 10;
  Gmin[1] = Gmax[1] =  0;
  solver.setInput(Gmin,"lbg");
  solver.setInput(Gmax,"ubg");

  // Solve the problem
  solver.solve();

  // Get the solution
  solver.getOutput(Usol,"x");
  cout << "optimal solution: " << Usol << endl;

  }
  
  return 0;

}

