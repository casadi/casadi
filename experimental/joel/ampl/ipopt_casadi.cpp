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
#include <casadi/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <nonlinear_programming/symbolic_nlp.hpp>

using namespace CasADi;

int main(int argc, char **argv){
  // Get the problem
  casadi_assert(argc==2);
  std::string problem = argv[1];
  
  // Create an NLP instance
  SymbolicNLP nlp;

  // Parse an NL-file
  nlp.parseNL(problem);
  
  // Objective function
  SXFunction ffcn(nlp.x,nlp.f);

  // NLP solver
  IpoptSolver nlp_solver;
  
  // Check if constrained
  if(nlp.g.size()>0){

    // Constraint function
    SXFunction gfcn(nlp.x,nlp.g);
    nlp_solver = IpoptSolver(ffcn,gfcn);
  } else {
    nlp_solver = IpoptSolver(ffcn);
  }
  
  // Set options
//   nlp_solver.setOption("max_iter",10);
//  nlp_solver.setOption("verbose",true);
//  nlp_solver.setOption("linear_solver","ma57");
  nlp_solver.setOption("generate_hessian",true);
//   nlp_solver.setOption("hessian_approximation","limited-memory");
//   nlp_solver.setOption("derivative_test","second-order");
  
  // Initialize NLP solver
  nlp_solver.init();
  
  // Pass the bounds and initial guess
  nlp_solver.setInput(nlp.x_lb,NLP_LBX);
  nlp_solver.setInput(nlp.x_ub,NLP_UBX);
  nlp_solver.setInput(nlp.g_lb,NLP_LBG);
  nlp_solver.setInput(nlp.g_ub,NLP_UBG);
  nlp_solver.setInput(nlp.x_init,NLP_X_INIT);
  
  // Solve NLP
  nlp_solver.solve();
  
  return 0;
}
