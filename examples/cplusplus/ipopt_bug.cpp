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
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/std_vector_tools.hpp>

using namespace CasADi;
using namespace std;

int main(){
  cout << "program started" << endl;

  SX x = SX::sym("x",1,1);

  // Create the NLP
  SXFunction nlp(nlpIn("x",x),nlpOut("f",x*x));

  // Allocate an NLP solver
  IpoptSolver solver(nlp);

  // initialize the solver
  solver.setOption("linear_solver","ma57");
  // solver.setOption("fixed_variable_treatment","make_constraint");  // THIS FIXES IT
  solver.init();

  // Bounds on u and initial condition
  vector<double> xmin(1), xmax(1), xinit(1);
  xmin[0] = 0;
  xmax[0] = 0;
  xinit[0] = 0;

  solver.setInput(xmin,"lbx");
  solver.setInput(xmax,"ubx");
  solver.setInput(xinit,"x0");

  // Solve the problem
  solver.solve();

  // Print the optimal cost
  double cost;
  solver.getOutput(cost,"f");
  cout << "optimal cost: " << cost << endl;

  // Print the optimal solution
  vector<double> xopt(1);
  solver.getOutput(xopt,"x");
  cout << "optimal x: " << xopt << endl;

  return 0;
}
