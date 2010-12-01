/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
#include "casadi/stl_vector_tools.hpp"
#include "ipopt_interface/ipopt_solver.hpp"
#include "casadi/expression_tools.hpp"

using namespace CasADi;
using namespace std;

// Const function
template<typename T>
T L(T p, T u){
  return p*p + 0.1*u*u;
}

// Optimal control
template<typename T>
T u_star(T lambda_v){
  return fmin(0.2, fmax(-0.2, -5*lambda_v));
}

int main(){
  try{
  
  // Dimensions
  int N = 1000;

  // Constants
  double T = 40; // horizon length

  // PART 1: CONSTRUCT THE INTEGRATOR
  // Initial position
  SX p0 = 1; // position
  SX v0 = 1; // speed
  SX lambda_p0("lambda_p0"); // position multiplier
  SX lambda_v0("lambda_v0"); // speed multiplier

  SX dt = T/N; // time step

  // Integrate over the intervals with Euler forward
  SX p=p0, v=v0, lambda_p=lambda_p0, lambda_v=lambda_v0;
  for(int j=0; j<N; ++j){
    // Calculate the control
    SX uj = u_star(lambda_v);

    // Take a step in the integration
    p += dt*v;
    v += dt*(-p+uj);
    lambda_p += dt*(2*p-lambda_v);
    lambda_v += dt*lambda_p;
  }

  // Terminal constraints
  SXMatrix g;
  g << lambda_p << lambda_v;

  // Degrees of freedom
  SXMatrix x;
  x << lambda_p0 << lambda_v0;

#if 0

  // Objective function
  SXMatrix f = trans(g)*g;

  // Create the NLP
  NLP nlp(f,x);

  // Bounds on x and initial condition
  vector<double> xmin(2), xmax(2), xsol(2);
  for(int i=0; i<2; ++i){
    xmin[i] = -1000;
    xmax[i] =  1000;
    xsol[i] = 0;
  }
  xsol[0] = 6.765;
  xsol[1] = 6.581;

  nlp.setVariableBounds(&xmin[0],&xmax[0]);
  nlp.setPrimal(&xsol[0]);

  // Allocate an NLP solver
  IpoptSolver solver(nlp);

  // Set options
  solver.setOption("exact_hessian",false);
  solver.setOption("abstol",1e-10);

  // initialize the solver
  solver.init();

  // Solve the problem
  solver.solve();

  // Get the optimal cost
  double fopt = nlp.eval_f();
  cout << "optimal cost: " << fopt << endl;

  // Get the optimal solution
  nlp.getPrimal(&xsol[0]);
  cout << "optimal solution: " << xsol << endl;
  
#else

  // Create a function
  SXFunction fcn(x,g);

  // current value for x
  vector<double> xval(2);
  xval[0] = 6.765;
  xval[1] = 6.581;

  // value for g
  vector<double> gval(2);

  // perturbed g
  vector<double> gper(2);

  // Perturbance
  double eps = 0.01;

  // Jacobian of g
  double jacg[2][2];

  // Step
  vector<double> dx(2);

  for(int k=0; k<1; ++k){

  // Evaluate g
  for(int i=0; i<fcn.input().data().size(); ++i)
    fcn.input().data()[i] = xval[i];

  fcn.evaluate();
  fcn.output().getDense(&gval[0]);

  
  // Perturb in the two directions
  for(int i=0; i<2; ++i){
    xval[i] += eps;

    for(int r=0; r<fcn.input().data().size(); ++r)
      fcn.input().data()[r] = xval[r];

    fcn.evaluate();
    fcn.output().getDense(&gper[0]);
    xval[i] -= eps;

    // Add the corresponding entries to the jacobian
    jacg[0][i] = (gper[0]-gval[0])/eps;
    jacg[1][i] = (gper[1]-gval[1])/eps;
  }
  
  // Calculate the determinant of the jacobian
  double detj = jacg[0][0]*jacg[1][1] - jacg[1][0]*jacg[0][1];
  
  // Calculate the inverse
  double invj[2][2];
  invj[0][0] = jacg[1][1]/detj;
  invj[1][0] = -jacg[1][0]/detj;
  invj[0][1] = -jacg[0][1]/detj;
  invj[1][1] = jacg[0][0]/detj;

  // Calculate the step size
  dx[0] = invj[0][0]*(-gval[0]) + invj[0][1]*(-gval[1]);
  dx[1] = invj[1][0]*(-gval[0]) + invj[1][1]*(-gval[1]);

  // Update x
  xval[0] += 0.1*dx[0];
  xval[1] += 0.1*dx[1];
  }

  // output
  cout << "xval" << xval << endl;
  cout << "dx" << dx << endl;
  cout << "gval " << gval << endl;
  cout << "jacg" << "{"<< jacg[0][0] << ", "<< jacg[0][1] << ", "<< jacg[1][0] << ", "<< jacg[1][1] << endl;
  
#endif  

  return 0;

  } catch (const char * str){
  cerr << str << endl;
  return 1;
}

  
}

