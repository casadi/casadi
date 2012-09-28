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
 *  Example program demonstrating sensitivity analysis with IPOPT from CasADi
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
  SXMatrix x = ssym("x",3);
  
  // Parameters
  SXMatrix p = ssym("p",2);
  
  // Objective
  SXMatrix f = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  
  // Constraints
  SXMatrix g = vertcat(
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
  vector<double> p_a  = {5.00,1.00};
  
  // Perturbed parameter values
  vector<double> p_b  = {4.50,1.00};

  // Initial guess and bounds for the optimization variables
  vector<double> x0  = {0.15, 0.15, 0.00, p_a[0], p_a[1]};
  vector<double> lbx = {0.00, 0.00, 0.00,   -inf,   -inf};
  vector<double> ubx = { inf,  inf,  inf,    inf,    inf};
  
  // Nonlinear bounds
  vector<double> lbg = {0.00, 0.00, p_a[0], p_a[1]};
  vector<double> ubg = {0.00, 0.00, p_a[0], p_a[1]};
    
  // Create NLP solver
  SXFunction ffcn(x,f);
  SXFunction gfcn(x,g);
  IpoptSolver solver(ffcn,gfcn);
  
  // Mark the parameters amongst the constraints (see sIPOPT documentation)
  Dictionary con_integer_md;
  con_integer_md["sens_init_constr"] = vector<int>{0,0,1,2};
  solver.setOption("con_integer_md",con_integer_md);
  
  // Mark the parameters amongst the variables (see sIPOPT documentation)
  Dictionary var_integer_md;
  var_integer_md["sens_state_1"] = vector<int>{0,0,0,1,2};
  solver.setOption("var_integer_md",var_integer_md);

  // Pass the perturbed values (see sIPOPT documentation)
  Dictionary var_numeric_md;
  var_numeric_md["sens_state_value_1"] = vector<double>{0,0,0,p_b[0],p_b[1]};
  solver.setOption("var_numeric_md",var_numeric_md);
  
  // Enable sensitivities
  solver.setOption("run_sens","yes");
  solver.setOption("n_sens_steps", 1);
  
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
  cout << "f_opt = " << solver.output(NLP_COST) << endl;
  cout << "x_opt = " << solver.output(NLP_X_OPT) << endl;
  
  return 0;
}

