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
#include <casadi/casadi.hpp>
#include <nonlinear_programming/sqp_method.hpp>
#include <interfaces/qpoases/qpoases_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <nonlinear_programming/lifted_sqp.hpp>

// Initialize the NLP at the optimal solution
bool inititialize_at_solution = false;

// Initialize the NLP at the optimal solution
bool inititialize_at_measurements = true;

// Initialize the NLP at the optimal solution
bool with_ipopt = false;

// Use Gauss-Newton
// bool gauss_newton = true;

// Lifted
bool lifted = true;

// Single-shooting
bool single_shooting = false;

// Expand MX->SX
// bool expand = true;

using namespace CasADi;
using namespace std;

int main(){

  // Physical parameters
  double g = 9.81; // gravity
  double poolwidth = 0.2;
  double sprad = 0.03;
  double spheight = 0.01;
  double endtime = 1.0;

  // Optimization parameters
  double drag_true = 2.0; // => u(0)
  double depth_true = 0.01; // => u(1)
  
  // Discretization

  // The largest dimensions which work with SX and IPOPT
  int numboxes = 3;
  int num_eulersteps = 20;
  int num_measurements = 20;

  // The largest dimensions which work with SX and exact Hessian
//   int numboxes = 20; // (but e.g. 10 fails with the error message "Initial QP could not be solved due to unboundedness"
//   int num_eulersteps = 10;
//   int num_measurements = 25;

  // The largest dimensions which work with SX and Gauss-Newton Hessian
  //   int numboxes = 20;
  //   int num_eulersteps = 10;
  //   int num_measurements = 50;

  // Plotting
  bool plot_progress = false;
  int numboxes_per_plot = 1;

  // Discretization
  int ntimesteps = num_eulersteps*num_measurements;
  double dt = endtime/ntimesteps;
  double dx = poolwidth/numboxes;
  double dy = poolwidth/numboxes;
  vector<double> x(numboxes), y(numboxes);
  for(int i=0; i<numboxes; ++i){
    x[i] = (i+0.5)*dx;
    y[i] = (i+0.5)*dy;
  }

  // Initial conditions
  DMatrix u0 = DMatrix::zeros(numboxes+1,numboxes  );
  DMatrix v0 = DMatrix::zeros(numboxes  ,numboxes+1);
  DMatrix h0 = DMatrix::zeros(numboxes  ,numboxes  ); 
  for(int i=0; i<numboxes; ++i){
    for(int j=0; j<numboxes; ++j){
      double spdist = sqrt(pow((x[i]-0.04),2.) + pow((y[j]-0.04),2.));
      if(spdist<sprad/3.0){
	h0.elem(i,j) = spheight * cos(3.0*M_PI*spdist/(2.0*sprad));
      }
    }
  }
  
  // Free parameters
  SX drag("b");
  SX depth("H");
  vector<SX> p(2); p[0]=drag; p[1]=depth;
  vector<double> p_true(2); p_true[0]=drag_true; p_true[1]=depth_true;
  
  // The state at a measurement
  SXMatrix uk = ssym("uk",numboxes+1, numboxes);
  SXMatrix vk = ssym("vk",numboxes  , numboxes+1);
  SXMatrix hk = ssym("hk",numboxes  , numboxes);
  
  // Take one step of the integrator
  SXMatrix u = uk;
  SXMatrix v = vk;
  SXMatrix h = hk;
  
  // Temporaries
  SX d1 = -dt*g/dx;
  SX d2 = dt*drag;
  
  // Update u
  for(int i=0; i<numboxes-1; ++i){
    for(int j=0; j<numboxes; ++j){
      u.elem(1+i,j) += d1*(h.elem(1+i,j)-h.elem(i,j))- d2*u.elem(1+i,j);
    }
  }
  
  // Update v
  d1 = -dt*g/dy;
  for(int i=0; i<numboxes; ++i){
    for(int j=0; j<numboxes-1; ++j){
      v.elem(i,j+1) += d1*(h.elem(i,j+1)-h.elem(i,j))- d2*v.elem(i,j+1);
    }
  }
    
  // Update h
  d1 = (-depth*dt)*(1.0/dx);
  d2 = (-depth*dt)*(1.0/dy);
  for(int i=0; i<numboxes; ++i){
    for(int j=0; j<numboxes; ++j){
      h.elem(i,j) += d1*(u.elem(1+i,j)-u.elem(i,j)) + d2*(v.elem(i,j+1)-v.elem(i,j));
    }
  }
  
  SXFunction f_step;
  {
    // Create an integrator function
    vector<SXMatrix> f_step_in(4);
    f_step_in[0] = p;
    f_step_in[1] = uk;
    f_step_in[2] = vk;
    f_step_in[3] = hk;
    vector<SXMatrix> f_step_out(3);
    f_step_out[0] = u;
    f_step_out[1] = v;
    f_step_out[2] = h;
    f_step = SXFunction(f_step_in,f_step_out);
//     f_step.setOption("live_variables",true);
    f_step.init();
    cout << "generated single step dynamics (" << f_step.getAlgorithmSize() << " nodes)" << endl;
  }

  // Integrate over one interval
  vector<MX> f_in(4);
  MX P = msym("P",2);
  MX Uk = msym("Uk",numboxes+1, numboxes);
  MX Vk = msym("Vk",numboxes  , numboxes+1);
  MX Hk = msym("Hk",numboxes  , numboxes);
  f_in[0] = P;
  f_in[1] = Uk;
  f_in[2] = Vk;
  f_in[3] = Hk;
  vector<MX> f_inter = f_in;
  vector<MX> f_out;
  for(int j=0; j<num_eulersteps; ++j){
    // Create a call node
    f_out = f_step.call(f_inter);
    
    // Save intermediate state
    f_inter[1] = f_out[0];
    f_inter[2] = f_out[1];
    f_inter[3] = f_out[2];
  }

  // Create an integrator function
  FX f = MXFunction(f_in,f_out);
  f.init();
  cout << "generated discrete dynamics, MX (" << shared_cast<MXFunction>(f).countNodes() << " nodes)" << endl;

  // Expand the discrete dynamics
  if(false){
    f = SXFunction(shared_cast<MXFunction>(f));
    f.init();
    cout << "generated discrete dynamics, SX (" << shared_cast<SXFunction>(f).getAlgorithmSize() << " nodes)" << endl;
  }

  // Measurements
  vector<DMatrix> H_meas;
  H_meas.reserve(num_measurements);

  // Simulate once to generate "measurements"
  f.setInput(p_true,0);
  f.setInput(u0,1);
  f.setInput(v0,2);
  f.setInput(h0,3);
  clock_t time1 = clock();
  for(int k=0; k<num_measurements; ++k){
    f.evaluate();
    const DMatrix& u = f.output(0);
    const DMatrix& v = f.output(1);
    const DMatrix& h = f.output(2);
    f.setInput(u,1);
    f.setInput(v,2);
    f.setInput(h,3);
    
    // Save a copy of h
    H_meas.push_back(h);
  }
  clock_t time2 = clock();
  double t_elapsed = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "measurements generated in " << t_elapsed << " seconds." << endl;

  // Initial guesses
  vector<double> drag_guess, depth_guess;
  drag_guess.push_back( 0.5); depth_guess.push_back(0.01);
  drag_guess.push_back( 5.0); depth_guess.push_back(0.01);
  drag_guess.push_back(15.0); depth_guess.push_back(0.01);
  drag_guess.push_back(30.0); depth_guess.push_back(0.01);
  drag_guess.push_back( 2.0); depth_guess.push_back(0.005);
  drag_guess.push_back( 2.0); depth_guess.push_back(0.02);
  drag_guess.push_back( 2.0); depth_guess.push_back(0.1);
  drag_guess.push_back( 0.2); depth_guess.push_back(0.001);
  drag_guess.push_back( 1.0); depth_guess.push_back(0.005);
  drag_guess.push_back( 4.0); depth_guess.push_back(0.02);
  drag_guess.push_back( 1.0); depth_guess.push_back(0.02);
  drag_guess.push_back(20.0); depth_guess.push_back(0.001);
  
  // Number of tests
  const int n_tests = drag_guess.size();

  // Number of iterations
  vector<int> iter_count_gn(n_tests,-1);
  vector<int> iter_count_eh(n_tests,-1);
  
  // Solution time
  vector<double> sol_time_gn(n_tests,-1);
  vector<double> sol_time_eh(n_tests,-1);
  
  // Estimated drag and depth
  vector<double> drag_est_gn(n_tests,-1);
  vector<double> depth_est_gn(n_tests,-1);
  vector<double> drag_est_eh(n_tests,-1);
  vector<double> depth_est_eh(n_tests,-1);
  
  for(int sol=0; sol<2; ++sol){
    bool gauss_newton = sol==0;
  
    // Number of NLP variables
    int nX = 2;
    if(!single_shooting) 
      nX += numboxes*numboxes*num_measurements;
    
    // Variables in the lifted NLP
    MX X = msym("X",nX);
    P = X[Slice(0,2)];
    int X_offset = 2;

    // Least-squares objective function
    MX F;
    
    // Constraint function
    MX G = MX::sparse(0,1);
    
    // Generate full-space NLP
    MX U = u0;
    MX V = v0;
    MX H = h0;
    int offset = 0;
    for(int k=0; k<num_measurements; ++k){
      // Take a step
      MX f_arg[4] = {P,U,V,H};
      vector<MX> f_res = f.call(vector<MX>(f_arg,f_arg+4));
      U = f_res[0];
      V = f_res[1];
      H = f_res[2];
      
      if(!single_shooting){
	// Lift the variable
	MX H_def = H;
	H = X[Slice(X_offset,X_offset+numboxes*numboxes)];
	H = reshape(H,h0.sparsity());
	X_offset += H.size();
	
	// Constraint function term
	G.append(flatten(H_def-H));
      }

      // Objective function term
      F.append(flatten(H-H_meas[k]));
    }

    // Number of lifted
    int nv = G.size();

    // Append variables
//     G.append(P);

	// Least-squares objective
    if(!gauss_newton){
      F = inner_prod(F,F)/2;
    }

    // Function which calculates the objective terms and constraints
    vector<MX> fg_out;
    fg_out.push_back(F);
    fg_out.push_back(G);
    MXFunction fg(X,fg_out);
    fg.init();
    cout << "Generated lifted NLP (" << fg.countNodes() << " nodes)" << endl;
    
    // Expand NLP
    SXFunction fg_expanded(fg);
    fg_expanded.init();
    cout << "expanded lifted NLP (" << fg_expanded.getAlgorithmSize() << " nodes)" << endl;
    
    // Formulate the NLP
    SXMatrix nlp_x = fg_expanded.inputSX();
    SXMatrix nlp_f = fg_expanded.outputSX(0);
    SXMatrix nlp_g = fg_expanded.outputSX(1);
    SXFunction ffcn(nlp_x,nlp_f);
    SXFunction gfcn(nlp_x,nlp_g);

    // Solve with IPOPT
    NLPSolver nlp_solver;
    if(with_ipopt){
      SXFunction hfcn;
      if(gauss_newton){
	// Scalar objective function
	ffcn = SXFunction(nlp_x,inner_prod(nlp_f,nlp_f)/2);
	
	// Gauss-Newton Hessian approximation for Ipopt
	SXMatrix JF = jacobian(nlp_f,nlp_x);
	hfcn = SXFunction(nlp_x,mul(trans(JF),JF));
	hfcn.setOption("live_variables",true);
      }
      nlp_solver = IpoptSolver(ffcn,gfcn,hfcn);
      //     if(!gauss_newton) nlp_solver.setOption("generate_hessian",true);
    } else {
      nlp_solver = LiftedSQP(ffcn,gfcn);
      if(gauss_newton) nlp_solver.setOption("gauss_newton",true);
      nlp_solver.setOption("qp_solver",Interfaces::QPOasesSolver::creator);
      Dictionary qp_solver_options;
      qp_solver_options["printLevel"] = "none";
      //qp_solver_options["verbose"] = true;
      nlp_solver.setOption("qp_solver_options",qp_solver_options);
      if(lifted) nlp_solver.setOption("num_lifted",nv);
      nlp_solver.setOption("maxiter",100);
      nlp_solver.setOption("toldx",1e-9);
      nlp_solver.setOption("verbose",true);
    }
    nlp_solver.init();
    
    // Run tests
    for(int test=0; test<n_tests; ++test){
      // Print progress
      cout << "test " << test << endl;
      try{
      
	// Initial guess for the parameters
	nlp_solver.input(NLP_X_INIT).setZero();
	if(inititialize_at_solution){
	  nlp_solver.setInput(p_true,NLP_X_INIT);
	} else {
	  nlp_solver.input(NLP_X_INIT)[0] = drag_guess[test];
	  nlp_solver.input(NLP_X_INIT)[1] = depth_guess[test];
	}

	// Initial guess for the heights
	if(!single_shooting && (inititialize_at_measurements || inititialize_at_solution)){
	  vector<double>::iterator it=nlp_solver.input(NLP_X_INIT).begin()+2;
	  for(int k=0; k<num_measurements; ++k){
	    copy(H_meas[k].begin(),H_meas[k].end(),it);
	    it += numboxes*numboxes;
	  }
	}
	nlp_solver.input(NLP_LBG).setZero();
	nlp_solver.input(NLP_UBG).setZero();
// 	nlp_solver.input(NLP_UBG)[nv] = numeric_limits<double>::infinity();
// 	nlp_solver.input(NLP_UBG)[nv+1] = numeric_limits<double>::infinity();
// 	nlp_solver.input(NLP_LBG)[nv] = -numeric_limits<double>::infinity();
// 	nlp_solver.input(NLP_LBG)[nv+1] = -numeric_limits<double>::infinity();
	
      //   nlp_solver.input(NLP_LBX).setAll(-1000.);
      //   nlp_solver.input(NLP_UBX).setAll( 1000.);
        nlp_solver.input(NLP_LBX)[0] = 0;
        nlp_solver.input(NLP_LBX)[1] = 0;

	time1 = clock();
	nlp_solver.solve();
	time2 = clock();
	
	// Solution statistics
	int& iter_count   = gauss_newton ? iter_count_gn[test] : iter_count_eh[test];
	double& sol_time  = gauss_newton ? sol_time_gn[test] : sol_time_eh[test];
	double& drag_est  = gauss_newton ? drag_est_gn[test] : drag_est_eh[test];
	double& depth_est = gauss_newton ? depth_est_gn[test] : depth_est_eh[test];
	
	iter_count = nlp_solver.getStat("iter_count");
	sol_time = double(time2-time1)/CLOCKS_PER_SEC;
	drag_est = nlp_solver.output(NLP_X_OPT).at(0);
	depth_est = nlp_solver.output(NLP_X_OPT).at(1);
	
	// Estimated drag
      } catch(exception& ex){
	cout << "Test " << test << " failed: " << ex.what() << endl;
      }
    }
  }

  // Tolerance 
  double tol=1e-3;
  
  cout << 
  setw(10) << "drag" <<  "  &" <<
  setw(10) << "depth" << "  &" << 
  setw(10) << "iter_gn" << "  &" << 
  setw(10) << "time_gn" << "  &" <<
  setw(10) << "iter_en" << "  &" << 
  setw(10) << "time_eh" << "  \\\\ \%" <<
  setw(10) << "edrag_gn" << 
  setw(10) << "edepth_gn" <<
  setw(10) << "edrag_eh" << 
  setw(10) << "edepth_eh" << endl;
  for(int test=0; test<n_tests; ++test){
    cout << setw(10) << drag_guess[test] << "  &";
    cout << setw(10) << depth_guess[test] << "  &";
    if(fabs(drag_est_gn[test]-drag_true) + fabs(depth_est_gn[test]-depth_true)<tol){
      cout << setw(10) << iter_count_gn[test] << "  &";
      cout << setw(10) << sol_time_gn[test] << "  &";
    } else {
      cout << setw(10) << "$\\infty$" << "  &";
      cout << setw(10) << "$\\infty$" << "  &";
    }
    if(fabs(drag_est_eh[test]-drag_true) + fabs(depth_est_eh[test]-depth_true)<tol){
      cout << setw(10) << iter_count_eh[test] << "  &";
      cout << setw(10) << sol_time_eh[test] << "  \\\\ \%";
    } else {
      cout << setw(10) << "$\\infty$" << "  &";
      cout << setw(10) << "$\\infty$" << "  \\\\ \%";
    }
    cout << setw(10) << drag_est_gn[test];
    cout << setw(10) << depth_est_gn[test];
    cout << setw(10) << drag_est_eh[test];
    cout << setw(10) << depth_est_eh[test] << endl;
  }
  
  return 0;
}
