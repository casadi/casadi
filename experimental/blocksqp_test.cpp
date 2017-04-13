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


/** \brief Writing a multiple shooting code from scratch

  This example demonstrates how to write a simple, yet powerful multiple shooting code from
  scratch using CasADi. For clarity, the code below uses a simple formulation with only states
  and controls, and no path constrants. It relies on CasADi machinery to keep track of sparsity,
  formulate ODE sensitivity equations and build up the Jacobian of the NLP constraint function.

  By extending the code below, it should be possible for a user to solve ever more complex
  problems. For example, one can easily make a multi-stage formulation by simply allocating
  another integrator instance and use the two integrators instances for different shooting nodes.
  By replacing the explicit CVodes integrator with the fully-implicit IDAS integrator, one can
  solve optimal control problems in differential-algebraic equations.

  \author Joel Andersson
  \date 2011-2012
*/

// CasADi core
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std;

int main(){
  auto x=SX::sym("x");
  auto y=SX::sym("y");
  SXDict nlp = {{"x", SX::vertcat({x,y})},
                {"f", sq(1-x)+100*sq(y-sq(x))},
                {"g", x+y}};
  Dict solver_options;
  auto solver = nlpsol("mysolver", "blocksqp", nlp, solver_options);
  DMDict solver_in;
  solver_in["x0"]=vector<double>{0,1};
  solver_in["lbx"]=vector<double>{-10,1.2};
  solver_in["ubx"]=vector<double>{10,2};
  solver_in["lbg"]=-10;
  solver_in["ubg"]=10;
  auto solver_out = solver(solver_in);
  
  /*
  // another simple test example
  auto x1=SX::sym("x1");
  auto x2=SX::sym("x2");
  auto x3=SX::sym("x3");
  auto x4=SX::sym("x4");
  auto x5=SX::sym("x5");
  SXDict nlp = {{"x", SX::vertcat({x1,x2,x3,x4,x5})},
                {"f", x1*sq(x2)+x3*sq(x4)+sq(x5)},
                {"g", sq(x1)+x2-x3+pow(x4,3)+x5}};
  Dict solver_options;
  solver_options["hess_update"] = 4;
  solver_options["max_conv_qp"] = 4;
  //solver_options["conv_strategy"] = 1;
  auto solver = nlpsol("mysolver", "blocksqp", nlp, solver_options);
  DMDict solver_in;
  solver_in["x0"]=vector<double>{0.5,1,1.5,2,1};
  solver_in["lbx"]=vector<double>{-10,-10,-10,-10,-10};
  solver_in["ubx"]=vector<double>{10,10,10,10,10};
  solver_in["lbg"]=10;
  solver_in["ubg"]=10;
  auto solver_out = solver(solver_in);
  cout << "opt = " << solver_out.at("x") << endl;
  */

  return 0;
}
