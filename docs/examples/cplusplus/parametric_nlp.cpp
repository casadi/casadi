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
 *  Example program demonstrating parametric NLPs in CasADi
 *  Note that there is currently no support for parametric sensitivities via this feature (although it would make a lot of sense).
 *  For parametric sensitivities, see the parametric_sensitivities.cpp example which calculates sensitivitied via the sIPOPT extension
 *  to IPOPT.
 * 
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
  
  // Infinity
  double inf = numeric_limits<double>::infinity();
  
  // Initial guess and bounds for the optimization variables
  double x0_[]  = {0.15, 0.15, 0.00};
  double lbx_[] = {0.00, 0.00, 0.00};
  double ubx_[] = { inf,  inf,  inf};
  vector<double> x0(x0_,x0_+3);
  vector<double> lbx(lbx_,lbx_+3);
  vector<double> ubx(ubx_,ubx_+3);
  
  // Nonlinear bounds
  double lbg_[] = {0.00, 0.00};
  double ubg_[] = {0.00, 0.00};
  vector<double> lbg(lbg_,lbg_+2);
  vector<double> ubg(ubg_,ubg_+2);
  
  // Original parameter values
  double p0_[]  = {5.00,1.00};
  vector<double> p0(p0_,p0_+2);

  // NLP
  SXFunction nlp(nlpIn("x",x,"p",p),nlpOut("f",f,"g",g));

  // Create NLP solver
  NlpSolver solver("ipopt", nlp);
  
  // Set options and initialize solver
  solver.init();
  
  // Solve NLP
  solver.setInput( x0, "x0");
  solver.setInput( p0, "p");
  solver.setInput(lbx, "lbx");
  solver.setInput(ubx, "ubx");
  solver.setInput(lbg, "lbg");
  solver.setInput(ubg, "ubg");
  solver.evaluate();
  
  // Print the solution
  cout << "-----" << endl;
  cout << "Optimal solution for p = " << str(solver.input("p")) << ":" << endl;
  cout << setw(30) << "Objective: " << str(solver.output("f")) << endl;
  cout << setw(30) << "Primal solution: " << str(solver.output("x")) << endl;
  cout << setw(30) << "Dual solution (x): " << str(solver.output("lam_x")) << endl;
  cout << setw(30) << "Dual solution (g): " << str(solver.output("lam_g")) << endl;
  
  // Change the parameter and resolve
  p0[0] = 4.5;
  solver.setInput( p0, "p");
  solver.evaluate();
  
  // Print the new solution
  cout << "-----" << endl;
  cout << "Optimal solution for p = " << str(solver.input("p")) << ":" << endl;
  cout << setw(30) << "Objective: " << str(solver.output("f")) << endl;
  cout << setw(30) << "Primal solution: " << str(solver.output("x")) << endl;
  cout << setw(30) << "Dual solution (x): " << str(solver.output("lam_x")) << endl;
  cout << setw(30) << "Dual solution (g): " << str(solver.output("lam_g")) << endl;
  
  return 0;
}

