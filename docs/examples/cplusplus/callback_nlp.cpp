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
  This example demonstrates how a plain C function can be embedded into an
  MX computational graph.
  
  Note that the syntax is not yet mature and will change in future releases.
 

 \author Joris Gillis
 \date 2015
*/

using namespace casadi;

// Plain C function to be embedded
double myFunction(double arg1, double arg2) {
  return (arg1-arg2)*(arg1-arg2);
}

// Wrapper Callback class
class MyFunction : public Callback2 {

  // The method that casadi will call to numerically evaluate the Wrapper
  virtual std::vector<DMatrix> operator()(const std::vector<DMatrix>& arg) {
    double a = arg[0].getValue();
    double b = arg[1].getValue();
    return std::vector<DMatrix>(1, myFunction(a,b));
  };
  
  // Our function has two inputs
  virtual int nIn() { return 2;}
};

int main(int argc, char **argv){

  MX x = MX::sym("x");
  MX y = MX::sym("y");
  
  // Create a pure CasADi Function out of the wrapper class
  MyFunction wrapper = MyFunction();
  Function fun = wrapper.create();
  
  // Call the pure CasADi function symbolically
  std::vector<MX> args = {2*x,cos(y)};
  MX z = fun(args)[0];
  
  MX vars = vertcat(x,y);

  // NLP
  MXFunction nlp("nlp", nlpIn("x",vars),nlpOut("f",z));

  // Set options
  Dict opts;
  
  // Note: we never specified derivatives for our C function.
  // CasADi falls back on finite differences.
  // This means we better not ask for exact Hessian,
  // and we cannot require a too small tolerance.
  opts["tol"] = 1e-5;
  opts["hessian_approximation"] = "limited-memory";

  // Allocate NLP solver and buffers
  NlpSolver nlp_solver("nlp_solver", "ipopt", nlp, opts);
  std::map<std::string, DMatrix> arg, res;
  
  res = nlp_solver(arg);
  
  std::cout << res["x"] << std::endl;
  return 0;
}
