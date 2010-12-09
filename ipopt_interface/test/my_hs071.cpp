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
#include "ipopt_interface/ipopt_solver.hpp"
#include "casadi/fx/sx_function.hpp"

using namespace CasADi;

int main(){

  try{

  // optimization variables
  SXMatrix x("x",4); // control

  // Objective function
  SXMatrix f = x[0]*x[3]*(x[0] + x[1] + x[2])  +  x[2];

  // Bounds on u
  std::vector<double> xmin(4), xmax(4);
  for(int i=0; i<4; ++i){
    xmin[i] = 1; 
    xmax[i] = 5;
  }

  // Constraints
  SXMatrix g(2,1);  std::vector<double> lb(2);   std::vector<double> ub(2);
  g[0] = x[0]*x[1]*x[2]*x[3];
  lb[0] = 25;  ub[0] = 2e19;
  
  g[1] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3];
  lb[1] = ub[1] = 40;
  
  // Starting point:
  std::vector<double> x0(4);
  x0[0]=1; x0[1]=5; x0[2]=5; x0[3]=1;

  // Create an NLP
  SXFunction ffcn(x,f);
  SXFunction gfcn(x,f);

  // Allocate an NLP solver
  IpoptSolver solver(ffcn,gfcn);

  // initialize the solver
  solver.init();

  // Pass the bounds
  solver.input(NLP_LBX).set(xmax);
  solver.input(NLP_UBX).set(xmax);
  solver.input(NLP_LBG).set(lb);
  solver.input(NLP_UBG).set(ub);

  // Pass initial guess
  solver.input(NLP_X_INIT).set(x0);
 
  // Solve the problem
  solver.solve();
  
  return 0;

  } catch (const char * str){
  std::cerr << str << std::endl;
  return 1;
}

  
}
