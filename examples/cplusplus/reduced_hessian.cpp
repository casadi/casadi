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
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>

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
  double x0_[]  = {25,0,0};
  double lbx_[] = {-inf, -inf, -inf};
  double ubx_[] = { inf,  inf,  inf};
  vector<double> x0(x0_,x0_+3);
  vector<double> lbx(lbx_,lbx_+3);
  vector<double> ubx(ubx_,ubx_+3);

  // Nonlinear bounds
  double lbg_[] = {0.00};
  double ubg_[] = {0.00};
  vector<double> lbg(lbg_,lbg_+1);
  vector<double> ubg(ubg_,ubg_+1);

  // Create NLP solver
  SXFunction nlp(nlpIn("x",x),nlpOut("f",f,"g",g));
  IpoptSolver solver(nlp);

  // Mark the parameters amongst the variables (see sIPOPT documentation)
  Dictionary var_integer_md;
  int red_hessian[] = {0,1,2};
  var_integer_md["red_hessian"] = std::vector<int>(red_hessian,red_hessian+3);
  solver.setOption("var_integer_md",var_integer_md);

  // Enable reduced hessian calculation
  solver.setOption("compute_red_hessian","yes");
  
  // Initialize solver
  solver.init();

  // Solve NLP
  solver.setInput( x0, "x0");
  solver.setInput(lbx, "lbx");
  solver.setInput(ubx, "ubx");
  solver.setInput(lbg, "lbg");
  solver.setInput(ubg, "ubg");
  solver.evaluate();

  // Print the solution
  cout << "----" << endl;
  cout << "Minimal cost " << solver.output("f") << endl;
  cout << "----" << endl;

  cout << "Solution" << endl;
  cout << "x = " << solver.output("x").data() << endl;
  cout << "----" << endl;
  
  // Get the reduced Hessian
  try{
    DMatrix red_hess = solver.getReducedHessian();
    cout << "Reduced Hessian = " << endl;
    red_hess.printDense();
  } catch(...){
    cout << "Support for retrieving the reduced Hessian not enabled." << endl;
  }
  
  return 0;
}

