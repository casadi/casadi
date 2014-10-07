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
#include <fstream>

using namespace casadi;
using namespace std;
/**
 *  Solve an NLP using codegen  
 *  Part 1: generation
 *  Joel Andersson, K.U. Leuven 2013
 */

void generateCode(Function fcn, const std::string& name, bool expand, std::ostream& makefile){
  cout << "Generating code for " << name << endl;

  // Convert to an SXFunction (may or may not improve efficiency)
  if(expand && is_a<MXFunction>(fcn)){
    fcn = SXFunction(shared_cast<MXFunction>(fcn));
    fcn.init();
  }

  // Generate C code
  fcn.generateCode(name + ".c");
  
  // Generate compilation instructions
  makefile << "add_library(" << name << " SHARED " << name << ".c)" << endl;
  makefile << "set_target_properties(" << name << " PROPERTIES PREFIX \"\")" << endl;
  makefile << "set_target_properties(" << name << " PROPERTIES SUFFIX \".casadi\")" << endl;
  makefile << endl;
}

int main(){
    
  /** Test problem 
   * 
   *    min x0^2 + x1^2
   *    s.t.    x0 + x1 - 10 = 0
   */

  // Optimization variables
  MX x = MX.sym("x",2);

  // Objective
  MX f = x[0]*x[0] + x[1]*x[1];

  // Constraints
  MX g = x[0]+x[1]-10;
    
  // Infinity
  double inf = numeric_limits<double>::infinity();

  // Convert MXFunction to SXFunction before code generation (may or may not improve efficiency)
  bool expand = true;

  // NLP function
  Function nlp = MXFunction(nlpIn("x",x),nlpOut("f",f,"g",g));
  nlp.init();

  // Gradient of the Lagrangian
  Function grad_f = nlp.gradient("x","f");
  grad_f.init();

  // Jacobian of the constraints
  Function jac_g = nlp.jacobian("x","g");
  jac_g.init();

  // Hessian of the lagrangian
  Function grad_lag = nlp.derivative(0,1);
  Function hess_lag = grad_lag.jacobian(NL_X,NL_NUM_OUT+NL_X,false,true);
  hess_lag.init();

  // Generate Makefile
  ofstream makefile;
  makefile.open("./CMakeLists.txt");
  makefile << "cmake_minimum_required(VERSION 2.6)" << endl;
  makefile << "project(nlp-codegen-autogen C)" << endl;

  // Codegen and compile
  generateCode(nlp,"nlp", expand, makefile);
  generateCode(grad_f,"grad_f", expand, makefile);
  generateCode(jac_g,"jac_g", expand, makefile);
  generateCode(hess_lag,"hess_lag", expand, makefile);

  // Finalize makefile
  makefile.close();

  return 0;
}

