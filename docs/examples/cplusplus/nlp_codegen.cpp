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
 *  Example program demonstrating NLP solution with Ipopt with callback functions as generated code
 *  Joel Andersson, 2013-2016
 */

int main(){
  /** Test problem
   *
   *    min x0^2 + x1^2
   *    s.t.    x0 + x1 - 10 = 0
   */

  // Optimization variables
  MX x = MX::sym("x", 2);

  // Objective
  MX f = x(0)*x(0) + x(1)*x(1);

  // Constraints
  MX g = x(0)+x(1)-10;

  // Create an NLP solver instance
  Function solver = nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}});

  // Generate C code for the NLP functions
  solver.generate_dependencies("nlp.c");

  // Just-in-time compilation?
  bool jit = false;
  if (jit) {
    // Create a new NLP solver instance using just-in-time compilation
    solver = nlpsol("solver", "ipopt", "nlp.c");
  } else {
    // Compile the c-code
    int flag = system("gcc -fPIC -shared -O3 nlp.c -o nlp.so");
    casadi_assert(flag==0, "Compilation failed");

    // Create a new NLP solver instance from the compiled code
    solver = nlpsol("solver", "ipopt", "nlp.so");
  }

  // Bounds and initial guess
  std::map<std::string, DM> arg, res;
  arg["lbx"] = -DM::inf();
  arg["ubx"] =  DM::inf();
  arg["lbg"] =  0;
  arg["ubg"] =  DM::inf();
  arg["x0"] = 0;

  // Solve the NLP
  res = solver(arg);

  // Print solution
  cout << "-----" << endl;
  cout << "objective at solution = " << res.at("f") << endl;
  cout << "primal solution = " << res.at("x") << endl;
  cout << "dual solution (x) = " << res.at("lam_x") << endl;
  cout << "dual solution (g) = " << res.at("lam_g") << endl;

  return 0;
}
