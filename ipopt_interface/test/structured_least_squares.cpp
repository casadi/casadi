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
#include <ipopt_interface/ipopt_solver.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>

using namespace CasADi;

// Philippe Dreesen's structured least squares problem from the numerical optimization course

int main(){
  try{

  // Hankel matrix
  double A_data[7][3] = {{-1,7,-1},{7,-1,10},{-1,10,4},{10,4,9},{4,9,1},{9,1,4},{1,4,3}};
  vector<double> A_vec(&A_data[0][0], &A_data[0][0]+7*3);
  SXMatrix A(A_vec,7,3);
  std::cout << "A = " << A << std::endl;

  // Approximation error in compact form
  SXMatrix e("e",9), v("v",3);

  // Tv matrix
  SXMatrix Tv(7,9);
  for(int i=0; i<7; ++i){
    Tv(i,i)   = v[0];
    Tv(i,i+1) = v[1];
    Tv(i,i+2) = v[2];
  }
  std::cout << "Tv = " << Tv << std::endl;

  // Objective function
  SXMatrix f = trans(e)*e/2;

  // Bounds on e
  SXMatrix emin = -10 * ones(e.size1(),e.size2());
  SXMatrix emax = 10 * ones(e.size1(),e.size2());

  // Bounds on v
  SXMatrix vmin = -10 * ones(v.size1(),v.size2());
  SXMatrix vmax = 10 * ones(v.size1(),v.size2());

  // Constraints
  SXMatrix g;
  g << trans(v)*v-1;
  g << A*v-Tv*e;
  std::cout << "g = " << g << std::endl;

  // Bounds on g
  SXMatrix gmin(g.size1(),g.size2()), gmax(g.size1(),g.size2());

  // All variables
  SXMatrix x, xmin, xmax;
  x << e; xmin << emin; xmax << emax;
  x << v; xmin << vmin; xmax << vmax;
  std::cout << "x = " << x << std::endl;
  std::cout << "xmin = " << xmin << std::endl;
  std::cout << "xmax = " << xmax << std::endl;

  // Initial guess
  int q = v.numel();
  std::vector<double> x0(x.numel());
  for(int i=0; i<e.numel(); ++i)
    x0[i] = 0.1;
  for(int i=e.numel(); i<x0.size(); ++i)
    x0[i] = 1/sqrt(q);

  std::cout << "x0 = " << x0 << std::endl;  

  // NLP
  SXFunction ffcn(x,f);
  SXFunction gfcn(x,g);

  // Bounds (cleanup)
  std::vector<double> xmin_val(xmin.numel());
  std::vector<double> xmax_val(xmax.numel());
  std::vector<double> gmin_val(gmin.numel());
  std::vector<double> gmax_val(gmax.numel());
  getValue(xmin,&xmin_val[0]);
  getValue(xmax,&xmax_val[0]);
  getValue(gmin,&gmin_val[0]);
  getValue(gmax,&gmax_val[0]);

  // Allocate an SQP solver
  IpoptSolver solver(ffcn,gfcn);

  // Set options
  solver.setOption("exact_hessian",true);
  solver.setOption("abstol",1e-10);

  // initialize the solver
  solver.init();

  // Set the starting point for the iteration
  solver.input(NLP_X_INIT).set(x0);
  solver.input(NLP_LBX).set(xmin_val);
  solver.input(NLP_UBX).set(xmax_val);
  solver.input(NLP_LBG).set(gmin_val);
  solver.input(NLP_UBG).set(gmax_val);
  
  // Solve the problem
  solver.solve();
  
  // Get the solution
  std::vector<double> xopt(x0.size());
  solver.output(NLP_X_OPT).get(xopt);
  std::cout << "optimal solution: " << xopt << std::endl;

  return 0;

  } catch (const char * str){
  std::cerr << str << std::endl;
  return 1;
}

  
}