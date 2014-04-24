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

#include <casadi/interfaces/liftopt/liftopt_solver.hpp>
#include <casadi/interfaces/liftopt/liftopt_internal.hpp>
#include <casadi/core/function/sx_function.hpp>
#include <casadi/core/mx/mx_tools.hpp>
#include <casadi/core/std_vector_tools.hpp>

using namespace casadi;
using namespace std;

int main(){
  const int timeSteps = 50;
  const double stepSize = 3.0 / timeSteps;

  // Declare some scalar variables
  SX t = SX::sym("t"), u = SX::sym("u");

  // Formulate the integrator
  vector<SX> F_in(2);
  F_in[0] = t;
  F_in[1] = u;
  SXFunction F(F_in, t + stepSize*(t*(t+1) + u));
  F.init();
  
  // Free variable
  MX U = MX::sym("U",timeSteps);

  // Build up a graph of calls to the function F - this graph will be lifted!
  MX trajectoryVal = 0.08;
  vector<MX> objResv(timeSteps*2);
  for(int k = 0; k < timeSteps; ++k){
    vector<MX> F_arg(2);
    F_arg[0] = trajectoryVal;
    F_arg[1] = U[k];
    trajectoryVal = F.call(F_arg)[0];
    objResv[k] = U[k];
    objResv[timeSteps+k] = trajectoryVal;
  }
  
  // Equality constraint
  MX eq = trajectoryVal;
  
  // Objective residual
  MX objRes = vertcat(objResv);
  
  // Objective value
  MX obj = inner_prod(objRes,objRes);

  // Lagrange multiplier
  MX lambda("lambda",1);
  
  // Lagrangian
  MX lag = obj - inner_prod(lambda,eq);

  // NLP cost and objective function argument
  vector<MX> fcn_in(LO_NUM_IN);
  fcn_in[LO_U] = U;
  fcn_in[LO_LAMBDA] = lambda;
  
  // NLP cost and objective function
  vector<MX> fcn_out(LO_NUM_OUT);
  fcn_out[LO_OBJRES] = objRes;
  fcn_out[LO_EQ] = eq;
  fcn_out[LO_INEQ] = MX(0,0); // FIXME: MX should be a 0-by-0 matrix by default, or?
  fcn_out[LO_OBJ] = obj;
  fcn_out[LO_LAGFCN] = lag;
  MXFunction fcn(fcn_in,fcn_out);
  
  // Create a liftopt solver instance
  LiftoptSolver solver(fcn);
  
  // Set some options
  solver.setOption("optimizer","sqp");
  solver.setOption("lifted",true);

  // Initialize the solver
  solver.init();
  
  // Pass inputs
  solver.input("x0").setAll(0);
  solver.input("lam_x0").setAll(0);
  solver.input("lbx").setAll(-1);
  solver.input("ubx").setAll(1);

  // Solve the problem
  solver.solve();
  
  return 0;
}
