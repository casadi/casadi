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

using namespace casadi;
using namespace std;

/**
Solve the Rosenbrock problem, formulated as the NLP:

minimize     x^2+tanh(y)^2
 x, y
subject to   cos(x+y)+0.5 == 0
             sin(x) + 0.5 <= 0

Joris Gillis, 2018
*/

int main(){

  // Declare variables
  SX x = SX::sym("x");
  SX y = SX::sym("y");

  // Formulate the NLP
  SX f = pow(x, 2) + pow(tanh(y), 2);
  SX g = vertcat(cos(x+y)+0.5, sin(x)+0.5);
  SXDict nlp = {{"x", SX::vertcat({x,y})},
                {"f", f},
                {"g", g}};

  Dict opts;
  opts["outer_iterations"] = 12;
  // Create an NLP solver
  Function solver = nlpsol("solver", "panoc", nlp, opts);

  // Solve the Rosenbrock problem
  DMDict arg;
  arg["x0"] = vector<double>{-0.5, -1.8};
  arg["lbg"] = vector<double>{0, -inf};
  arg["ubg"] = vector<double>{0, 0};
  DMDict res = solver(arg);

  //  Print solution
  cout << "Optimal cost:                     " << res.at("f") << endl;
  cout << "Primal solution:                  " << res.at("x") << endl;
  cout << "Dual solution (simple bounds):    " << res.at("lam_x") << endl;
  cout << "Dual solution (nonlinear bounds): " << res.at("lam_g") << endl;

  cout << "Iteration count:" << solver.stats()["iter_count"] << endl;

  return 0;
}
