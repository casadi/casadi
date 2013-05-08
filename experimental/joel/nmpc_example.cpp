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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSefcn.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <symbolic/casadi.hpp>
#include <interfaces/qpoases/qpoases_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <nonlinear_programming/nlp_qp_solver.hpp>
#include <nonlinear_programming/scpgen.hpp>

#include <iomanip>
#include <ctime>
#include <cstdlib>

using namespace CasADi;
using namespace std;

int main(){

  // Horizon length
  double tf = 3.0;
  
  // Number of subintervals
  int n = 30;
  
  // Time step
  MX dt = tf/n;
  
  // Parameter (should be treated as such)
  double x0 = 0.02;
  
  // Control
  MX u = msym("u",n);
  vector<double> lbu(n, -1.0);
  vector<double> ubu(n,  1.0);

  // Objective function terms
  vector<MX> ff;
    
  // Add control regularization
  ff.push_back(u);

  // Constraints
  vector<MX> gg;
  vector<double> lbg, ubg;

  // Lifting modes
  enum LiftMode{UNLIFTED, AUT_INIT, ZERO_INIT};
  LiftMode mode = ZERO_INIT;

  // Perform lifted single-shooting
  MX x = x0;
  for(int k=0; k<n; ++k){
    // Integrate
    x = x + dt*(x*(x+1) + u[k]);

    // Lift the state
    switch(mode){
    case AUT_INIT: x.lift(x); break;
    case ZERO_INIT: x.lift(0.); break;
    case UNLIFTED: break;
    }

    // Objective function terms
    ff.push_back(x);

    // State bounds
    gg.push_back(x);
    if(k==n-1){
      lbg.push_back( 0.0);
      ubg.push_back( 0.0);
    } else {
      lbg.push_back(-1.0);
      ubg.push_back( 1.0);
    }
  }

  // Gather least square terms and constraints
  MX f = vertcat(ff);
  MX g = vertcat(gg);
    
  // Use Gauss-Newton?
  bool gauss_newton = true;
  if(!gauss_newton){
    f = inner_prod(f,f)/2;
  }

  // Form the NLP
  MXFunction nlp(nlIn("x",u),nlOut("f",f,"g",g));
  SCPgen solver(nlp);

  //solver.setOption("verbose",true);
  solver.setOption("regularize",false);
  solver.setOption("codegen",false);
  solver.setOption("maxiter_ls",1);
  solver.setOption("maxiter",100);
  if(gauss_newton){
    solver.setOption("hessian_approximation","gauss-newton");
  }
  
  // Print the variables
  solver.setOption("print_x",range(0,n,5));

  Dictionary qp_solver_options;
  if(false){
    solver.setOption("qp_solver",NLPQPSolver::creator);
    qp_solver_options["nlp_solver"] = IpoptSolver::creator;
    Dictionary nlp_solver_options;
    nlp_solver_options["tol"] = 1e-12;
    nlp_solver_options["print_level"] = 0;
    nlp_solver_options["print_time"] = false;
    qp_solver_options["nlp_solver_options"] = nlp_solver_options;
      
  } else {
    solver.setOption("qp_solver",QPOasesSolver::creator);
    qp_solver_options["printLevel"] = "none";
  }
  solver.setOption("qp_solver_options",qp_solver_options);

  solver.init();
    
  // Pass bounds and solve
  solver.setInput(lbu,"lbx");
  solver.setInput(ubu,"ubx");
  solver.setInput(lbg,"lbg");
  solver.setInput(ubg,"ubg");
  solver.solve();

  cout << "u_opt = " << solver.output(NLP_SOLVER_X).data() << endl;

    
  return 0;
}

