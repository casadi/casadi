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
#include <casadi/casadi.hpp>

using namespace casadi;

/**
Solve the Rosenbrock problem, formulated as the NLP:

minimize     x^2 + 100*z^2
subject to   z+(1-x)^2-y == 0

Joel Andersson, 2015-2016
*/

int main(){

  // Declare variables
  SX x = SX::sym("x");
  SX y = SX::sym("y");
  SX z = SX::sym("z");

  // Formulate the NLP
  SX f = pow(x,2) + 100*pow(z,2);
  SX g = z + pow(1-x, 2) - y;
  SXDict nlp = {{"x", SX::vertcat({x,y,z})},
                {"f", f},
                {"g", g}};

  // Create an NLP solver
  Function solver = nlpsol("solver", "ipopt", nlp);

  // Solve the Rosenbrock problem
  DMDict arg;
  arg["x0"] = std::vector<double>{2.5,3.0,0.75};
  arg["lbg"] = arg["ubg"] = 0;
  DMDict res = solver(arg);

  //  Print solution
  std::cout << "Optimal cost:                     " << double(res.at("f")) << std::endl;
  std::cout << "Primal solution:                  " << std::vector<double>(res.at("x")) << std::endl;
  std::cout << "Dual solution (simple bounds):    " << std::vector<double>(res.at("lam_x")) << std::endl;
  std::cout << "Dual solution (nonlinear bounds): " << std::vector<double>(res.at("lam_g")) << std::endl;

  return 0;
}
