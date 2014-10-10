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
#include <casadi/core/misc/symbolic_nlp.hpp>

/**
 * This example demonstrates how NL-files, which can be generated
 * by AMPl or Pyomo, can be imported in CasADi and solved using
 * e.g. the interface to AMPL
 
 \author Joel Andersson, Vyacheslav Kungurtsev
 \date 2013
*/
 
 
using namespace casadi;
 
int main(int argc, char **argv){
  // Get the problem
  std::string problem = (argc==2) ? argv[1] : "../docs/examples/nl_files/hs107.nl";

  // Parse an NL-file
  SymbolicNLP nl;
  nl.parseNL(problem);

  // NLP
  SXFunction nlp(nlpIn("x",nl.x),nlpOut("f",nl.f,"g",nl.g));
 
  // Allocate NLP solver
  NlpSolver nlp_solver("sqpmethod", nlp);

  // Set options
  // nlp_solver.setOption("max_iter",10);
  // nlp_solver.setOption("verbose",true);
  // nlp_solver.setOption("linear_solver","ma57");
  nlp_solver.setOption("hessian_approximation","exact");
  // nlp_solver.setOption("derivative_test","second-order");

  // Specify QP solver
  nlp_solver.setOption("qp_solver","nlp");
  Dictionary qp_solver_options;
  qp_solver_options["nlp_solver"] = "ipopt"; 
  Dictionary nlp_solver_options;
  nlp_solver_options["print_level"] = 0;
  nlp_solver_options["print_time"] = 0;
  qp_solver_options["nlp_solver_options"] = nlp_solver_options;
  nlp_solver.setOption("qp_solver_options",qp_solver_options);
 
  // Initialize NLP solver
  nlp_solver.init();

  // Pass the bounds and initial guess
  nlp_solver.setInput(nl.x_lb,"lbx");
  nlp_solver.setInput(nl.x_ub,"ubx");
  nlp_solver.setInput(nl.g_lb,"lbg");
  nlp_solver.setInput(nl.g_ub,"ubg");
  nlp_solver.setInput(nl.x_init,"x0");

  // Solve NLP
  nlp_solver.evaluate();
  return 0;
}
