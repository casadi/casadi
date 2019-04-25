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
  arg["x0"] = vector<double>{2.5,3.0,0.75};
  arg["lbg"] = arg["ubg"] = 0;
  DMDict res = solver(arg);

  //  Print solution
  cout << "Optimal cost:                     " << double(res.at("f")) << endl;
  cout << "Primal solution:                  " << vector<double>(res.at("x")) << endl;
  cout << "Dual solution (simple bounds):    " << vector<double>(res.at("lam_x")) << endl;
  cout << "Dual solution (nonlinear bounds): " << vector<double>(res.at("lam_g")) << endl;

  return 0;
}
