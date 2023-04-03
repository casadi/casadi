/*
 *    MIT No Attribution
 *
 *    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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

/** Solves the following NLP with block structured Hessian:
    minimize      x^2 - 0.5*y^2
    subject to    x - y = 0

    Joel Andersson, 2016
*/

int main(){

  // Declare variables
  SX x = SX::sym("x");
  SX y = SX::sym("y");

  // Formulate the NLP
  SX f = pow(x,2) - 0.5*pow(y,2);
  SX g = x - y;
  SXDict nlp = {{"x", SX::vertcat({x,y})},
                {"f", f},
                {"g", g}};

  // Create an NLP solver
  Dict opts;
  opts["opttol"] = 1.0e-12;
  opts["nlinfeastol"] = 1.0e-12;
  opts["globalization"] = 0;
  opts["hess_update"] = 0;
  opts["hess_scaling"] = 0;
  opts["fallback_scaling"] = 0;
  opts["hess_lim_mem"] = 0;
  opts["max_consec_skipped_updates"] = 200;
  opts["block_hess"] = 0;
  opts["which_second_derv"] = 0;
  opts["schur"] = false;
  Function solver = nlpsol("solver", "blocksqp", nlp, opts);

  // Solve the Rosenbrock problem
  DMDict arg;
  arg["x0"] = std::vector<double>{10, 10};
  arg["lbg"] = arg["ubg"] = 0;
  DMDict res = solver(arg);

  //  Print solution
  std::cout << "Optimal cost:                     " << double(res.at("f")) << std::endl;
  std::cout << "Primal solution:                  " << std::vector<double>(res.at("x")) << std::endl;
  std::cout << "Dual solution (simple bounds):    " << std::vector<double>(res.at("lam_x")) << std::endl;
  std::cout << "Dual solution (nonlinear bounds): " << std::vector<double>(res.at("lam_g")) << std::endl;

  return 0;
}
