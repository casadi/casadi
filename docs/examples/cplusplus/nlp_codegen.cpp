/*
 *    MIT No Attribution
 *
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <casadi/casadi.hpp>

using namespace casadi;
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
  arg["ubg"] =  0;
  arg["x0"] = 0;

  // Solve the NLP
  res = solver(arg);

  // Print solution
  std::cout << "-----" << std::endl;
  std::cout << "objective at solution = " << res.at("f") << std::endl;
  std::cout << "primal solution = " << res.at("x") << std::endl;
  std::cout << "dual solution (x) = " << res.at("lam_x") << std::endl;
  std::cout << "dual solution (g) = " << res.at("lam_g") << std::endl;

  return 0;
}
