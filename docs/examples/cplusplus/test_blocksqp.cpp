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
  arg["x0"] = vector<double>{10, 10};
  arg["lbg"] = arg["ubg"] = 0;
  DMDict res = solver(arg);

  //  Print solution
  cout << "Optimal cost:                     " << double(res.at("f")) << endl;
  cout << "Primal solution:                  " << vector<double>(res.at("x")) << endl;
  cout << "Dual solution (simple bounds):    " << vector<double>(res.at("lam_x")) << endl;
  cout << "Dual solution (nonlinear bounds): " << vector<double>(res.at("lam_g")) << endl;

  return 0;
}
