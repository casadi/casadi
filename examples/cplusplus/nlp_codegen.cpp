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

void generateCodeAndCompile(FX fcn, const std::string& name, bool expand){
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

  // Objective scaling factor
  MX sigma = msym("sigma");
  
  // Lagrange multipliers
  MX lam = msym("lam",g.sparsity());

  // Form the lagrangian
  MX lag = sigma*f + inner_prod(lam,g);

  // Convert MXFunction to SXFunction before code generation (may or may not improve efficiency)
  bool expand = true;

  // Objective function
  MXFunction ffcn(x,f); 
  ffcn.init();
  generateCodeAndCompile(ffcn,"ffcn", expand);

  // Gradient of the Lagrangian
  FX grad_ffcn = ffcn.gradient();
  grad_ffcn.init();
  generateCodeAndCompile(grad_ffcn,"grad_ffcn", expand);

  // Constraint function
  MXFunction gfcn(x,g);
  gfcn.setOption("numeric_jacobian",false); // NOTE!
  gfcn.init();
  generateCodeAndCompile(gfcn,"gfcn", expand);

  // Jacobian of the constraints
  FX jac_gfcn = gfcn.jacobian();
  jac_gfcn.init();
  generateCodeAndCompile(jac_gfcn,"jac_gfcn", expand);

  // Lagrange function
  vector<MX> lfcn_in;
  lfcn_in.push_back(x);
  lfcn_in.push_back(lam);
  lfcn_in.push_back(sigma);
  MXFunction lfcn(lfcn_in,lag);
  lfcn.setOption("numeric_hessian",false); // NOTE!
  lfcn.init();

  // Hessian of the lagrangian
  FX hess_lfcn = lfcn.hessian();
  hess_lfcn.init();
  generateCodeAndCompile(hess_lfcn,"hess_lfcn", expand);

  // Load the generated functions into CasADi
  ExternalFunction ffcn_e("./ffcn.so");
  ExternalFunction grad_ffcn_e("./grad_ffcn.so");
  ExternalFunction gfcn_e("./gfcn.so");
  ExternalFunction jac_gfcn_e("./jac_gfcn.so");
  ExternalFunction hess_lfcn_e("./hess_lfcn.so");

  // Create an NLP solver passing derivative information
  IpoptSolver solver(ffcn_e, gfcn_e, hess_lfcn_e, jac_gfcn_e, grad_ffcn_e);
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

