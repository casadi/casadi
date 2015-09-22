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
    
  /** Test problem (Ganesh & Biegler, A reduced Hessian strategy for sensitivity analysis of optimal flowsheets, AIChE 33, 1987, pp. 282-296)
   * 
   *    min     x1^2 + x2^2 + x3^2
   *    s.t.    6*x1 + 3&x2 + 2*x3 - pi = 0
   *            p2*x1 + x2 - x3 - 1 = 0
   *            x1, x2, x3 >= 0
   * 
   */
  
  // Optimization variables
  SX x = SX::sym("x",3);
  
  // Parameters
  SX p = SX::sym("p",2);
  
  // Objective
  SX f = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  
  // Constraints
  SX g = vertcat(
       6*x[0] + 3*x[1] + 2*x[2] - p[0],
    p[1]*x[0] +   x[1] -   x[2] -    1
  );
  
  // Augment the parameters to the list of variables
  x = vertcat(x,p);
  
  // Fix the parameters by additional equations
  g = vertcat(g,p);
  
  // Infinity
  double inf = numeric_limits<double>::infinity();
  
  // Original parameter values
  vector<double> p_a = {5.00, 1.00};
  
  // Perturbed parameter values
  vector<double> p_b = {4.50, 1.00};

  // Initial guess and bounds for the optimization variables
  vector<double> x0  = {0.15, 0.15, 0.00, p_a[0], p_a[1]};
  vector<double> lbx = {0.00, 0.00, 0.00,   -inf,   -inf};
  vector<double> ubx = { inf,  inf,  inf,    inf,    inf};
  
  // Nonlinear bounds
  vector<double> lbg = {0.00, 0.00, p_a[0], p_a[1]};
  vector<double> ubg = {0.00, 0.00, p_a[0], p_a[1]};
    
  // Create NLP
  SXFunction nlp("nlp", nlpIn("x", x), nlpOut("f", f, "g", g));

  // NLP solver options
  Dict opts;
  
  // Mark the parameters amongst the constraints (see sIPOPT documentation)
  Dict con_integer_md;
  con_integer_md["sens_init_constr"] = vector<int>{0, 0, 1, 2};
  opts["con_integer_md"] = con_integer_md;
  
  // Mark the parameters amongst the variables (see sIPOPT documentation)
  Dict var_integer_md;
  var_integer_md["sens_state_1"] = vector<int>{0, 0, 0, 1, 2};
  opts["var_integer_md"] = var_integer_md;

  // Pass the perturbed values (see sIPOPT documentation)
  Dict var_numeric_md;
  var_numeric_md["sens_state_value_1"] = vector<double>{0,0,0,p_b[0],p_b[1]};
  opts["var_numeric_md"] = var_numeric_md;
  
  // Enable sensitivities
  opts["run_sens"] = "yes";
  opts["n_sens_steps"] = 1;
  
  // Create NLP solver and buffers
  NlpSolver solver("solver", "ipopt", nlp, opts);
  std::map<std::string, DMatrix> arg, res;
  
  // Solve NLP
  arg["lbx"] = lbx;
  arg["ubx"] = ubx;
  arg["lbg"] = lbg;
  arg["ubg"] = ubg;
  arg["x0"] = x0;
  res = solver(arg);
  
  // Print the solution
  cout << "----" << endl;
  cout << "Minimal cost " << res.at("f") << endl;
  cout << "----" << endl;

  cout << "Nominal solution" << endl;
  cout << "x = " << res.at("x") << endl;
  cout << "----" << endl;
  
  cout << "perturbed solution" << endl;
  var_numeric_md = solver.getStat("var_numeric_md");
  cout << "x = " << var_numeric_md["sens_sol_state_1"] << endl;
  cout << "----" << endl;
  
  cout << "Dual bound multipliers" << endl;
  cout << "z_L = " << var_numeric_md["sens_sol_state_1_z_L"] << endl;
  cout << "z_U = " << var_numeric_md["sens_sol_state_1_z_U"] << endl;
  cout << "----" << endl;
  
  cout << "Constraint multipliers" << endl;
  Dict con_numeric_md = solver.getStat("con_numeric_md");
  cout << "lambda = " << con_numeric_md["sens_sol_state_1"] << endl;
  cout << "----" << endl;
  
  return 0;
}

