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
#include <fstream>
#include <ctime>
#include <iomanip>
#include <casadi/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <casadi/stl_vector_tools.hpp>

using namespace CasADi;
using namespace std;
/**
 *  Example program demonstrating sensitivity analysis with sIPOPT from CasADi
 *  Joel Andersson, K.U. Leuven 2012
 */

int main(){
    
  /** Test problem (from sIPOPT example collection)
   * 
   *    min  (x1-1)^2 +(x2-2)^2 + (x3-3)^2
   *    s.t.    x1+2*x2+3*x3 = 0
   * 
   */
  
  // Optimization variables
  SXMatrix x = ssym("x",3);
  
  // Objective
  SXMatrix f = pow(x[0]-1,2) + pow(x[1]-2,2) + pow(x[2]-3,2);
  
  // Constraint
  SXMatrix g = x[0]+2*x[1]+3*x[2];
  
  // Infinity
  double inf = numeric_limits<double>::infinity();

  // Initial guess and bounds for the optimization variables
  vector<double> x0  = {25,0,0};
  vector<double> lbx = {-inf, -inf, -inf};
  vector<double> ubx = { inf,  inf,  inf};

  // Nonlinear bounds
  vector<double> lbg = {0.00};
  vector<double> ubg = {0.00};

  // Create NLP solver
  SXFunction ffcn(x,f);
  SXFunction gfcn(x,g);
  IpoptSolver solver(ffcn,gfcn);

  // Mark the parameters amongst the variables (see sIPOPT documentation)
  Dictionary var_integer_md;
  var_integer_md["red_hessian"] = std::vector<int>{0,1,2};
  solver.setOption("var_integer_md",var_integer_md);

  // Enable reduced hessian calculation
  solver.setOption("compute_red_hessian","yes");
  
  // Initialize solver
  solver.init();

  // Solve NLP
  solver.setInput( x0, NLP_X_INIT);
  solver.setInput(lbx, NLP_LBX);
  solver.setInput(ubx, NLP_UBX);
  solver.setInput(lbg, NLP_LBG);
  solver.setInput(ubg, NLP_UBG);
  solver.evaluate();

  // Print the solution
  cout << "----" << endl;
  cout << "Minimal cost " << solver.output(NLP_COST) << endl;
  cout << "----" << endl;

  cout << "Solution" << endl;
  cout << "x = " << solver.output(NLP_X_OPT).data() << endl;
  cout << "----" << endl;
  
  return 0;
}

