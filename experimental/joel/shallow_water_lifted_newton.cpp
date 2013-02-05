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
#include <nonlinear_programming/lifted_sqp.hpp>

#include <iomanip>
#include <ctime>

using namespace CasADi;
using namespace std;

class Tester{
public:
  // Constructor
  Tester(int n_boxes, int n_euler, int n_meas) : n_boxes_(n_boxes), n_euler_(n_euler), n_meas_(n_meas){}

  // Perform the modelling
  void model();

  // Simulate to generae measurements
  void simulate(double drag_true, double depth_true);

  // Transscribe as an NLP
  void transcribe(bool single_shooting);

  // Perform the parameter estimation
  void optimize(double drag_guess, double depth_guess, int& iter_count, double& sol_time, double& drag_est, double& depth_est);

  // Dimensions
  int n_boxes_;
  int n_euler_;
  int n_meas_;

  // Initial conditions
  DMatrix u0_;
  DMatrix v0_;
  DMatrix h0_;

  // Discrete time dynamics
  FX f_;
  
  // Generated measurements
  vector<DMatrix> H_meas_;

  // NLP solver
  bool single_shooting_;
  NLPSolver nlp_solver_;    

};

void Tester::model(){
  // Physical parameters
  double g = 9.81; // gravity
  double poolwidth = 0.2;
  double sprad = 0.03;
  double spheight = 0.01;
  double endtime = 1.0;
    
  // Discretization
  int ntimesteps = n_euler_*n_meas_;
  double dt = endtime/ntimesteps;
  double dx = poolwidth/n_boxes_;
  double dy = poolwidth/n_boxes_;
  vector<double> x(n_boxes_), y(n_boxes_);
  for(int i=0; i<n_boxes_; ++i){
    x[i] = (i+0.5)*dx;
    y[i] = (i+0.5)*dy;
  }
  
  // Initial conditions
  u0_ = DMatrix::zeros(n_boxes_+1,n_boxes_  );
  v0_ = DMatrix::zeros(n_boxes_  ,n_boxes_+1);
  h0_ = DMatrix::zeros(n_boxes_  ,n_boxes_  ); 
  for(int i=0; i<n_boxes_; ++i){
    for(int j=0; j<n_boxes_; ++j){
      double spdist = sqrt(pow((x[i]-0.04),2.) + pow((y[j]-0.04),2.));
      if(spdist<sprad/3.0){
	h0_.elem(i,j) = spheight * cos(3.0*M_PI*spdist/(2.0*sprad));
      }
    }
  }
  
  // Free parameters
  SX drag("b");
  SX depth("H");
  vector<SX> p(2); p[0]=drag; p[1]=depth;
  
  // The state at a measurement
  SXMatrix uk = ssym("uk",n_boxes_+1, n_boxes_);
  SXMatrix vk = ssym("vk",n_boxes_  , n_boxes_+1);
  SXMatrix hk = ssym("hk",n_boxes_  , n_boxes_);
  
  // Take one step of the integrator
  SXMatrix u = uk;
  SXMatrix v = vk;
  SXMatrix h = hk;
  
  // Temporaries
  SX d1 = -dt*g/dx;
  SX d2 = dt*drag;
  
  // Update u
  for(int i=0; i<n_boxes_-1; ++i){
    for(int j=0; j<n_boxes_; ++j){
      u.elem(1+i,j) += d1*(h.elem(1+i,j)-h.elem(i,j))- d2*u.elem(1+i,j);
    }
  }
  
  // Update v
  d1 = -dt*g/dy;
  for(int i=0; i<n_boxes_; ++i){
    for(int j=0; j<n_boxes_-1; ++j){
      v.elem(i,j+1) += d1*(h.elem(i,j+1)-h.elem(i,j))- d2*v.elem(i,j+1);
    }
  }
  
  // Update h
  d1 = (-depth*dt)*(1.0/dx);
  d2 = (-depth*dt)*(1.0/dy);
  for(int i=0; i<n_boxes_; ++i){
    for(int j=0; j<n_boxes_; ++j){
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
  MX Uk = msym("Uk",n_boxes_+1, n_boxes_);
  MX Vk = msym("Vk",n_boxes_  , n_boxes_+1);
  MX Hk = msym("Hk",n_boxes_  , n_boxes_);
  f_in[0] = P;
  f_in[1] = Uk;
  f_in[2] = Vk;
  f_in[3] = Hk;
  vector<MX> f_inter = f_in;
  vector<MX> f_out;
  for(int j=0; j<n_euler_; ++j){
    // Create a call node
    f_out = f_step.call(f_inter);
    
    // Save intermediate state
    f_inter[1] = f_out[0];
    f_inter[2] = f_out[1];
    f_inter[3] = f_out[2];
  }
  
  // Create an integrator function
  f_ = MXFunction(f_in,f_out);
  f_.init();
  cout << "generated discrete dynamics, MX (" << shared_cast<MXFunction>(f_).countNodes() << " nodes)" << endl;
  
  // Expand the discrete dynamics
  if(false){
    f_ = SXFunction(shared_cast<MXFunction>(f_));
    f_.init();
    cout << "generated discrete dynamics, SX (" << shared_cast<SXFunction>(f_).getAlgorithmSize() << " nodes)" << endl;
  }
}

void Tester::simulate(double drag_true, double depth_true){
  
  // Measurements
  H_meas_.reserve(n_meas_);
  
  // Simulate once to generate "measurements"
  vector<double> p_true(2); p_true[0]=drag_true; p_true[1]=depth_true;
  f_.setInput(p_true,0);
  f_.setInput(u0_,1);
  f_.setInput(v0_,2);
  f_.setInput(h0_,3);
  clock_t time1 = clock();
  for(int k=0; k<n_meas_; ++k){
    f_.evaluate();
    const DMatrix& u = f_.output(0);
    const DMatrix& v = f_.output(1);
    const DMatrix& h = f_.output(2);
    f_.setInput(u,1);
    f_.setInput(v,2);
    f_.setInput(h,3);
    
    // Save a copy of h
    H_meas_.push_back(h);
  }
  clock_t time2 = clock();
  double t_elapsed = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "measurements generated in " << t_elapsed << " seconds." << endl;
}


void Tester::transcribe(bool single_shooting){
  single_shooting_ = single_shooting;
  
  // Number of NLP variables
  int nX = 2;
  if(!single_shooting) 
    nX += n_boxes_*n_boxes_*n_meas_;
  
  // Variables in the lifted NLP
  MX X = msym("X",nX);
  MX P = X[Slice(0,2)];
  int X_offset = 2;
  
  // Least-squares objective function
  MX F;
  
  // Constraint function
  MX G = MX::sparse(0,1);
  
  // Generate full-space NLP
  MX U = u0_;
  MX V = v0_;
  MX H = h0_;
  int offset = 0;
  for(int k=0; k<n_meas_; ++k){
    // Take a step
    MX f_arg[4] = {P,U,V,H};
    vector<MX> f_res = f_.call(vector<MX>(f_arg,f_arg+4));
    U = f_res[0];
    V = f_res[1];
    H = f_res[2];
    
    if(!single_shooting){
      // Lift the variable
      MX H_def = H;
      H = X[Slice(X_offset,X_offset+n_boxes_*n_boxes_)];
      H = reshape(H,h0_.sparsity());
      X_offset += H.size();
      
      // Constraint function term
      G.append(flatten(H_def-H));
    }
    
    // Objective function term
    F.append(flatten(H-H_meas_[k]));
  }
  
  // Number of lifted
  int nv = G.size();
  
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
  SXMatrix nlp_x = fg_expanded.inputExpr(0);
  SXMatrix nlp_f = fg_expanded.outputExpr(0);
  SXMatrix nlp_g = fg_expanded.outputExpr(1);
  SXFunction ffcn(nlp_x,nlp_f);
  SXFunction gfcn(nlp_x,nlp_g);
  
  
  
  nlp_solver_ = LiftedSQP(ffcn,gfcn);
  nlp_solver_.setOption("gauss_newton",true);
  nlp_solver_.setOption("qp_solver",QPOasesSolver::creator);
  Dictionary qp_solver_options;
  qp_solver_options["printLevel"] = "none";
  //qp_solver_options["verbose"] = true;
  nlp_solver_.setOption("qp_solver_options",qp_solver_options);
  nlp_solver_.setOption("num_lifted",nv);
  nlp_solver_.setOption("maxiter",100);
  nlp_solver_.setOption("toldx",1e-9);
  nlp_solver_.setOption("verbose",true);
  nlp_solver_.init();

}

void Tester::optimize(double drag_guess, double depth_guess, int& iter_count, double& sol_time, double& drag_est, double& depth_est){
  // Initial guess for the parameters
  nlp_solver_.input(NLP_X_INIT).setZero();
  nlp_solver_.input(NLP_X_INIT)[0] = drag_guess;
  nlp_solver_.input(NLP_X_INIT)[1] = depth_guess;

  // Initial guess for the heights
  if(!single_shooting_){
    vector<double>::iterator it=nlp_solver_.input(NLP_X_INIT).begin()+2;
    for(int k=0; k<n_meas_; ++k){
      copy(H_meas_[k].begin(),H_meas_[k].end(),it);
      it += n_boxes_*n_boxes_;
    }
  }
  
  nlp_solver_.input(NLP_LBG).setZero();
  nlp_solver_.input(NLP_UBG).setZero();
  // 	nlp_solver_.input(NLP_UBG)[nv] = numeric_limits<double>::infinity();
  // 	nlp_solver_.input(NLP_UBG)[nv+1] = numeric_limits<double>::infinity();
  // 	nlp_solver_.input(NLP_LBG)[nv] = -numeric_limits<double>::infinity();
  // 	nlp_solver_.input(NLP_LBG)[nv+1] = -numeric_limits<double>::infinity();
	
  //   nlp_solver_.input(NLP_LBX).setAll(-1000.);
  //   nlp_solver_.input(NLP_UBX).setAll( 1000.);
  nlp_solver_.input(NLP_LBX)[0] = 0;
  nlp_solver_.input(NLP_LBX)[1] = 0;

  clock_t time1 = clock();
  nlp_solver_.solve();
  clock_t time2 = clock();
  
  // Solution statistics  
  iter_count = nlp_solver_.getStat("iter_count");
  sol_time = double(time2-time1)/CLOCKS_PER_SEC;
  drag_est = nlp_solver_.output(NLP_X_OPT).at(0);
  depth_est = nlp_solver_.output(NLP_X_OPT).at(1);
}

int main(){

  double drag_true = 2.0; // => u(0)
  double depth_true = 0.01; // => u(1)
 

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
  
  // Create a tester object
  // Tester t(3,20,20); // The largest dimensions which work with SX and IPOPT
  Tester t(15,10,10); // The largest dimensions which work with SX and exact Hessian
  // Tester t(20,10,50); // The largest dimensions which work with SX and Gauss-Newton Hessian
    
  // Perform the modelling
  t.model();

  // Optimization parameters
  t.simulate(drag_true, depth_true);
  
  // For both single and multiple shooting
  for(int sol=0; sol<2; ++sol){

    // Transcribe as an NLP
    bool single_shooting = sol==0;
    t.transcribe(single_shooting);
    
    // Run tests
    for(int test=0; test<n_tests; ++test){
      // Print progress
      cout << "test " << test << endl;
      try{
	t.optimize(drag_guess[test],depth_guess[test],
		   sol==0 ? iter_count_gn[test] : iter_count_eh[test],
		   sol==0 ? sol_time_gn[test] : sol_time_eh[test],
		   sol==0 ? drag_est_gn[test] : drag_est_eh[test],
		   sol==0 ? depth_est_gn[test] : depth_est_eh[test]);

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
  setw(10) << "iter_ss" << "  &" << 
  setw(10) << "time_ss" << "  &" <<
  setw(10) << "iter_ms" << "  &" << 
  setw(10) << "time_ms" << "  \\\\ \%" <<
  setw(10) << "edrag_ss" << 
  setw(10) << "edepth_ss" <<
  setw(10) << "edrag_ms" << 
  setw(10) << "edepth_ms" << endl;
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
