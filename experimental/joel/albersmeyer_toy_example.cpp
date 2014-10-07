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

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <core/casadi.hpp>
#include <nonlinear_programming/sqp_method.hpp>
#include <interfaces/qpoases/qpoases_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <core/std_vector_tools.hpp>
#include <nonlinear_programming/lifted_sqp.hpp>

using namespace casadi;
using namespace std;

int main(){
  cout << "program started" << endl;
  
  // Initial values to test
  const int ntests = 11;
  double x0_test[ntests] = {0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30};

  // Solvers
  enum Solvers{IPOPT, LIFTED_SQP, FULLSPACE_SQP, OLD_SQP_METHOD, NUM_SOLVERS};
  string solver_name[] = {"IPOPT","LIFTED_SQP","FULLSPACE_SQP","OLD_SQP_METHOD"};

  // Iteration counts
  int iter_count[ntests][NUM_SOLVERS];
  
  // Lifting?
  bool lifting = true;
  
  // Gauss-Newton
  bool gauss_newton = false;
  
  // Tests for different initial values for x
  for(int test=0; test<ntests; ++test){
    double x0 = x0_test[test];
    cout << "x0 = " << x0 << endl;

    // Bounds on the state
    double x_init = 0;
    double x_min = -1;
    double x_max =  1;
    double xf_min = 0;
    double xf_max = 0;

    // End time
    double T = 3;

    // Dimensions
    int nk = 30;  // Number of control segments

    // Control
    vector<SX> u = ssym("u",nk).data(); // control

    // Time step
    SX dT = T/nk;

    // Intermediate variables with initial values and bounds
    vector<SX> v, v_eq;
    vector<double> v_init, v_min, v_max;

    // Objective terms
    vector<SX> F;

    // Get an expression for the state that the final time
    SX x = x0;
    for(int k=0; k<nk; ++k){
      // Get new value for X
      x = x + dT*(x*(x+1)+u[k]);
      
      // Append terms to objective function
      F.push_back(u[k]);
      F.push_back(x);
      
      if(lifting){
	// Lift x
	SX x_def = x;
	v_init.push_back(x_init);
	v_min.push_back(x_min);
	v_max.push_back(x_max);
	
	// Allocate intermediate variable
	stringstream ss;
	ss << "v_" << k;
	x = SX(ss.str());
	v.push_back(x);
	v_eq.push_back(x_def-x);
      }
    }
      
    // Objective
    SX f = 0;
    for(vector<SX>::const_iterator it=F.begin(); it!=F.end(); ++it){
      f += *it**it;
    }
    
    // Terminal constraints
    vector<SX> g(1,x);

    // Bounds on g
    vector<double> g_min(1,xf_min);
    vector<double> g_max(1,xf_max);

    // Bound and initial condition for the compact NLP
    vector<double> u_min(nk,-1);
    vector<double> u_max(nk, 1);
    vector<double> u_init(nk, 0);

    // Add intermediate variables to get the full-space NLP
    int nv = v.size(); // number of lifted variables
    u.insert(u.end(),v.begin(),v.end());
    u_min.insert(u_min.end(),v_min.begin(),v_min.end());
    u_max.insert(u_max.end(),v_max.begin(),v_max.end());
    u_init.insert(u_init.end(),v_init.begin(),v_init.end());
    g.insert(g.begin(),v_eq.begin(),v_eq.end());
    vector<double> v_bnd(nv,0);
    g_min.insert(g_min.begin(),v_bnd.begin(),v_bnd.end());
    g_max.insert(g_max.begin(),v_bnd.begin(),v_bnd.end());
    
    // Formulate the full-space NLP
    SXFunction ffcn;
    if(gauss_newton){
      ffcn = SXFunction(u,F);
    } else {
      ffcn = SXFunction(u,f);
    }
    SXFunction gfcn(u,g);

    Dictionary qp_solver_options;
    qp_solver_options["printLevel"] = "none";
    
    // Solve using multiple NLP solvers
    for(int solver=0; solver<NUM_SOLVERS; ++solver){
      cout << "Testing " << solver_name[solver] << endl;
      iter_count[test][solver] = 9999;
	
      // Get the nlp solver and NLP solver options
      NlpSolver nlp_solver;
      switch(solver){
	case IPOPT:
	  if(gauss_newton) continue; // not supported
	  nlp_solver = IpoptSolver(ffcn,gfcn);
	  if(gauss_newton){
	    nlp_solver.setOption("gauss_newton",true);
	  } else {
	    nlp_solver.setOption("generate_hessian",true);
	  }
	  nlp_solver.setOption("tol",1e-9);
	  break;
	case LIFTED_SQP:
	  if(!lifting && x0>=0.10) continue;;
	  nlp_solver = LiftedSQP(ffcn,gfcn);
	  if(gauss_newton)
	    nlp_solver.setOption("gauss_newton",true);
	  nlp_solver.setOption("qp_solver",QPOasesSolver::creator);
	  nlp_solver.setOption("qp_solver_options",qp_solver_options);
	  nlp_solver.setOption("num_lifted",nv);
	  nlp_solver.setOption("toldx",1e-9);
// 	  nlp_solver.setOption("verbose",true);
	  break;
	case FULLSPACE_SQP:
	  if(!lifting && x0>=0.10) continue;;
	  nlp_solver = LiftedSQP(ffcn,gfcn);
	  if(gauss_newton)
	    nlp_solver.setOption("gauss_newton",true);
	  nlp_solver.setOption("qp_solver",QPOasesSolver::creator);
	  nlp_solver.setOption("qp_solver_options",qp_solver_options);
	  nlp_solver.setOption("num_lifted",0);
	  nlp_solver.setOption("toldx",1e-9);
// 	  nlp_solver.setOption("verbose",true);
	  break;
	case OLD_SQP_METHOD:
	  if(gauss_newton) continue; // not supported
	  if(!lifting && x0>=0.07) continue;;
	  nlp_solver = SQPMethod(ffcn,gfcn);
	  nlp_solver.setOption("qp_solver",QPOasesSolver::creator);
	  nlp_solver.setOption("qp_solver_options",qp_solver_options);
	  if(gauss_newton){
	    nlp_solver.setOption("gauss_newton",true);
	  } else {
	    nlp_solver.setOption("generate_hessian",true);
	  }
      }
      
      // initialize the solver
      nlp_solver.init();

      // Initial guess and bounds
      nlp_solver.setInput(u_min,"lbx");
      nlp_solver.setInput(u_max,"ubx");
      nlp_solver.setInput(u_init,"x0");
      nlp_solver.setInput(g_min,"lbg");
      nlp_solver.setInput(g_max,"ubg");

      // Solve the problem
      nlp_solver.solve();
      
      // Get number of iterations
      iter_count[test][solver] = nlp_solver.getStat("iter_count");
      
      // Print the optimal solution
  //     cout << "optimal cost:    " << nlp_solver.output(NLP_SOLVER_F).toScalar() << endl;
  //     cout << "optimal control: " << nlp_solver.output(NLP_SOLVER_X) << endl;
  //     cout << "multipliers (u): " << nlp_solver.output(NLP_SOLVER_LAM_X) << endl;
  //     cout << "multipliers (gb): " << nlp_solver.output(NLP_SOLVER_LAM_G) << endl;
    }
  }
  
  // Make a table of the results
  cout << "---------------------------------------" << endl;
  for(int solver=0; solver<NUM_SOLVERS; ++solver){
    cout << setw(20) << solver_name[solver];
  }
  cout << endl;
  for(int test=0; test<ntests; ++test){
    for(int solver=0; solver<NUM_SOLVERS; ++solver){
      cout << setw(20) << iter_count[test][solver];
    }
    cout << endl;
  }
  return 0;
}




