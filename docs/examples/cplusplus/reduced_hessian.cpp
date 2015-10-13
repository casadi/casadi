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
#include <fstream>
#include <ctime>
#include <iomanip>
#include <casadi/casadi.hpp>

using namespace casadi;
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
  SX x = SX::sym("x",3);
  
  // Objective
  SX f = pow(x[0]-1,2) + pow(x[1]-2,2) + pow(x[2]-3,2);
  
  // Constraint
  SX g = x[0]+2*x[1]+3*x[2];
  
  // Infinity
  double inf = numeric_limits<double>::infinity();

  // Initial guess and bounds for the optimization variables
  vector<double> x0  = {25,0,0};
  vector<double> lbx = {-inf, -inf, -inf};
  vector<double> ubx = { inf,  inf,  inf};

  // Nonlinear bounds
  vector<double> lbg = {0.00};
  vector<double> ubg = {0.00};

  // Create NLP
  Function nlp = SX::fun("nlp", nlpIn("x", x), nlpOut("f", f, "g", g));

  // NLP solver options
  Dict opts;

  // Mark the parameters amongst the variables (see sIPOPT documentation)
  Dict var_integer_md;
  var_integer_md["red_hessian"] = std::vector<int>{0,1,2};
  opts["var_integer_md"] = var_integer_md;

  // Enable reduced hessian calculation
  opts["compute_red_hessian"] = "yes";
  
  // Create NLP solver and buffers
  NlpSolver solver("solver", "ipopt", nlp, opts);
  std::map<std::string, DMatrix> arg, res;

  // Solve NLP
  arg["x0"] = x0;
  arg["lbx"] = lbx;
  arg["ubx"] = ubx;
  arg["lbg"] = lbg;
  arg["ubg"] = ubg;
  res = solver(arg);

  // Print the solution
  cout << "----" << endl;
  cout << "Minimal cost " << res.at("f") << endl;
  cout << "----" << endl;

  cout << "Solution" << endl;
  cout << "x = " << res.at("x") << endl;
  cout << "----" << endl;
  
  // Get the reduced Hessian
  try{
    DMatrix red_hess = solver.getReducedHessian();
    cout << "Reduced Hessian = " << red_hess << endl;
  } catch(...){
    cout << "Support for retrieving the reduced Hessian not enabled." << endl;
  }
  
  return 0;
}

