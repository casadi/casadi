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

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>

using namespace CasADi;
using namespace std;
/**
 *  Example program demonstrating NLP solution with Ipopt with callback functions as generated code
 *  Joel Andersson, K.U. Leuven 2013
 */

FX generateCodeAndCompile(FX fcn, const std::string& name, bool expand){
  cout << "Generating code for " << name << endl;

  // Convert to an SXFunction (may or may not improve efficiency)
  if(expand && is_a<MXFunction>(fcn)){
    fcn = SXFunction(shared_cast<MXFunction>(fcn));
    fcn.init();
  }

  // Generate C code
  fcn.generateCode(name + ".c");

  // Compilation command
  string compile_command = "gcc -fPIC -shared -O3 " + name + ".c -o " + name + ".so";

  // Compile the c-code
  int flag = system(compile_command.c_str());
  casadi_assert_message(flag==0, "Compilation failed");

  // Load the generated function for evaluation
  ExternalFunction fcn_e("./" + name + ".so");
  return fcn_e;
}

int main(){
    
  /** Test problem 
   * 
   *    min x0^2 + x1^2
   *    s.t.    x0 + x1 - 10 = 0
   */

  // Optimization variables
  MX x = msym("x",2);

  // Objective
  MX f = x[0]*x[0] + x[1]*x[1];

  // Constraints
  MX g = x[0]+x[1]-10;
    
  // Infinity
  double inf = numeric_limits<double>::infinity();

  // Convert MXFunction to SXFunction before code generation (may or may not improve efficiency)
  bool expand = true;

  // NLP function
  FX nlp = MXFunction(nlpIn("x",x),nlpOut("f",f,"g",g));
  nlp.setOption("numeric_jacobian",false);
  nlp.init();

  // Gradient of the Lagrangian
  FX grad_f = nlp.gradient("x","f");
  grad_f.init();

  // Jacobian of the constraints
  FX jac_g = nlp.jacobian("x","g");
  jac_g.init();

  // Hessian of the lagrangian
  FX grad_lag = nlp.derivative(0,1);
  grad_lag.setOption("numeric_jacobian",false);
  FX hess_lag = grad_lag.jacobian(NLP_X,NLP_NUM_OUT+NLP_X,false,true);
  hess_lag.init();

  // Codegen and compile
  nlp = generateCodeAndCompile(nlp,"nlp", expand);
  grad_f = generateCodeAndCompile(grad_f,"grad_f", expand);
  jac_g = generateCodeAndCompile(jac_g,"jac_g", expand);
  hess_lag = generateCodeAndCompile(hess_lag,"hess_lag", expand);

  // Create an NLP solver passing derivative information
  IpoptSolver solver(nlp);
  solver.setOption("grad_f",grad_f);
  solver.setOption("jac_g",jac_g);
  solver.setOption("hess_lag",hess_lag);
  solver.init();

  // Set constraint bounds
  solver.setInput(0.,"lbg");

  // Solve the NLP
  solver.evaluate();

  // Print solution
  cout << "-----" << endl;
  cout << "objective at solution = " << solver.output("f") << endl;
  cout << "primal solution = " << solver.output("x") << endl;
  cout << "dual solution (x) = " << solver.output("lam_x") << endl;
  cout << "dual solution (g) = " << solver.output("lam_g") << endl;
  
  return 0;
}

