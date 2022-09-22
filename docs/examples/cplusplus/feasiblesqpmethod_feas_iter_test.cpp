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

#include <casadi/casadi.hpp>

/**
 * This example demonstrates how NL-files, which can be generated
 * by AMPl or Pyomo, can be imported in CasADi and solved using
 * e.g. the interface to AMPL

 \author Joel Andersson, Vyacheslav Kungurtsev
 \date 2013
*/

/*
DESCRIPTION: The problem is a simple QP. Any SQP solver should just need one
iteration! Depending on the size of the trust-region. It could take more 
iteration. But if the trust-region is big enough such that the solution is
included. The problem should just take one iteration!
*/
using namespace casadi;

int main(int argc, char **argv){
  
  // Set options
  Dict opts;
  opts["expand"] = true;
  // opts["max_iter"] = 10)
  // opts["verbose"] = true;
  // opts["linear_solver"] = "ma57";
  opts["hessian_approximation"] = "exact";
  opts["max_inner_iter"] = 5;
  // opts["derivative_test"] = "second-order";

  // Specify QP solver
  opts["qpsol"]  = "nlpsol";
  opts["qpsol_options.nlpsol"] = "ipopt";
  opts["qpsol_options.error_on_fail"] = false;
  opts["qpsol_options.nlpsol_options.ipopt.print_level"] = 0;
  opts["qpsol_options.nlpsol_options.ipopt.sb"] = "yes";
  opts["qpsol_options.nlpsol_options.print_time"] = 0;

  SX x = SX::sym("x");
  SX y = SX::sym("y");
  SX f = pow(x, 2) + pow(y-2, 2);;
  SX g = pow(x, 2) + pow(y, 2);
  SXDict nlp = {{"x", SX::vertcat({x,y})}, {"f", f}, {"g", g}};

  // Create an NLP solver
  // Function solver = nlpsol("solver", "sqpmethod", nlp, opts);
  Function solver = nlpsol("solver", "feasiblesqpmethod", nlp, opts);

  // Solve the problem
  DMDict arg;
  arg["x0"] = std::vector<double>{2.0, 0.0};//std::vector<double>{2.5,3.0,0.75};
  arg["lbg"] = 4;
  arg["ubg"] = 4;
  DMDict res = solver(arg);

  //  Print solution
  uout() << "Optimal cost:                     " << double(res.at("f")) << std::endl;
  uout() << "Primal solution:                  " << std::vector<double>(res.at("x")) << std::endl;
  return 0;
}
