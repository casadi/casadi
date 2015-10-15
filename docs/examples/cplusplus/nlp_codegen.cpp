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
 *  Joel Andersson, K.U. Leuven 2013
 */

Function generateCodeAndCompile(Function fcn, const std::string& name, bool expand){
  cout << "Generating code for " << name << endl;

  // Convert to sx (may or may not improve efficiency)
  if(expand && fcn.is_a("mxfunction")) {
    fcn = SX::fun(fcn.name(), fcn);
  }

  // Generate C code
  fcn.generate(name);

  // Compilation command
  string compile_command = "gcc -fPIC -shared -O3 " + name + ".c -o " + name + ".so";

  // Compile the c-code
  int flag = system(compile_command.c_str());
  casadi_assert_message(flag==0, "Compilation failed");

  // Load the generated function for evaluation
  return Function::external(name);
}

int main(){
  /** Test problem 
   * 
   *    min x0^2 + x1^2
   *    s.t.    x0 + x1 - 10 = 0
   */

  // Optimization variables
  MX x = MX::sym("x",2);

  // Objective
  MX f = x[0]*x[0] + x[1]*x[1];

  // Constraints
  MX g = x[0]+x[1]-10;
    
  // Convert mxfunction to sxfunction before code generation (may or may not improve efficiency)
  bool expand = true;

  // NLP function
  Function nlp = MX::fun("nlp", nlpIn("x", x),nlpOut("f", f, "g", g));

  // Gradient of the objective
  Function grad_f = nlp.gradient("x", "f");

  // Jacobian of the constraints
  Function jac_g = nlp.jacobian("x", "g");

  // Hessian of the lagrangian
  Function grad_lag = nlp.derivative(0, 1);
  Function hess_lag = grad_lag.jacobian(NL_X, NL_NUM_OUT+NL_X, false, true);

  // Codegen and compile
  nlp = generateCodeAndCompile(nlp,"nlp", expand);
  grad_f = generateCodeAndCompile(grad_f,"grad_f", expand);
  jac_g = generateCodeAndCompile(jac_g,"jac_g", expand);
  hess_lag = generateCodeAndCompile(hess_lag,"hess_lag", expand);

  // Create an NLP solver passing derivative information
  NlpSolver solver("solver", "ipopt", nlp,
                   Dict{{"grad_f", grad_f}, {"jac_g", jac_g}, {"hess_lag",hess_lag}});

  // Bounds and initial guess
  std::map<std::string, DMatrix> arg, res;
  arg["lbx"] = -DMatrix::inf();
  arg["ubx"] =  DMatrix::inf();
  arg["lbg"] =  0;
  arg["ubg"] =  DMatrix::inf();
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

