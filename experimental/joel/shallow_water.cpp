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

#include <iomanip>
#include <ctime>
#include <cstdlib>

using namespace CasADi;
using namespace std;

class Tester{
public:
  // Constructor
  Tester(int n_boxes, int n_euler, int n_finite_elements, int n_meas) : n_boxes_(n_boxes), n_euler_(n_euler), n_finite_elements_(n_finite_elements), n_meas_(n_meas){}

  // Perform the modelling
  void model();

  // Simulate to generae measurements
  void simulate(double drag_true, double depth_true);

  // Transscribe as an NLP
  void transcribe(bool single_shooting, bool gauss_newton);

  // Prepare NLP
  void prepare(bool codegen, bool ipopt_as_qp_solver, bool regularization, double reg_threshold);
 
  // Solve NLP
  void solve(int& iter_count);

  // Perform the parameter estimation
  void optimize(double drag_guess, double depth_guess, int& iter_count, double& sol_time, double& drag_est, double& depth_est);

  // Codegen function
  void dynamicCompilation(MXFunction& f, FX& f_gen, std::string fname, std::string fdescr);

  // Calculate the L1-norm of the primal infeasibility
  double primalInfeasibility();

  // Calculate the L1-norm of the dual infeasibility
  double dualInfeasibility();

  // Evaluate the residual function
  void eval_res();

  // Form the condensed QP
  void eval_qpf();

  // Regularize the condensed QP
  void regularize();

  // Solve the QP to get the (full) step
  void solve_qp();

  // Perform the line-search to take the step
  void line_search(int& ls_iter, bool& ls_success);

  // Evaluate the step expansion
  void eval_exp();

  // Print iteration header
  void printIteration(std::ostream &stream);
  
  // Print iteration
  void printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, 
                      double reg, int ls_trials, bool ls_success);

  // Dimensions
  int n_boxes_;
  int n_euler_;
  int n_finite_elements_;
  int n_meas_;

  // Initial conditions
  DMatrix u0_;
  DMatrix v0_;
  DMatrix h0_;

  // Discrete time dynamics
  FX f_;
  
  // Generated measurements
  vector<DMatrix> H_meas_;

  // Options
  bool single_shooting_;
  bool verbose_;
  bool gauss_newton_;
  bool regularization_;
  double reg_threshold_;

  // NLP
  SXFunction fg_sx_;
  MXFunction fg_mx_;

  /// QP solver for the subproblems
  QPSolver qp_solver_;

  /// maximum number of sqp iterations
  int maxiter_; 

  /// stopping criterion for the stepsize
  double toldx_;
  
  /// stopping criterion for the lagrangian gradient
  double tolgl_;
  
  /// Outputs of the linearization function
  enum QPFOut{QPF_G,QPF_H,QPF_B,QPF_A,QPF_NUM_OUT};
  
  /// Generate initial guess for lifted variables
  FX vinit_fcn_;

  /// Residual function
  FX res_fcn_;
 
  /// Quadratic approximation
  FX qp_fcn_;

  /// Step expansion
  FX exp_fcn_;
  
  /// Dimensions
  int nu_, ng_, ngL_;

  // Objective value
  double obj_k_;

  // Simple and nonlinear bounds
  vector<double> lbu_, ubu_, g_, lbg_, ubg_, gL_;

  /// Multipliers for the nonlinear bounds
  vector<double> lambda_g_, dlambda_g_;

  int res_lam_g_;
  int res_obj_, res_gl_, res_g_;

  int z_lam_g_;
  int z_obj_, z_gl_, z_g_;

  int qpf_lam_g_;

  int exp_du_, exp_dlam_g_, exp_lam_g_;
  int exp_osens_;

  struct Var{
    int n;

    int res_var, res_lam;
    int res_d, res_lam_d;

    int z_var, z_lam;
    int z_def, z_defL;

    int qpf_var, qpf_lam, qpf_res, qpf_resL;

    int exp_var, exp_lam;
    int exp_def, exp_defL;

    vector<double> step, init, opt, lam, dlam;
    vector<double> res, resL;
    
  };
  vector<Var> x_;  

  // Best merit function encountered so far
  double meritmax_;
  
  // Penalty parameter of merit function
  double sigma_;

  // 1-norm of last primal step
  double pr_step_;
  
  // 1-norm of last dual step
  double du_step_;

  // Regularization
  double reg_;

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
  bool any_point_in_domain = false;
  for(int i=0; i<n_boxes_; ++i){
    for(int j=0; j<n_boxes_; ++j){
      double spdist = sqrt(pow((x[i]-0.04),2.) + pow((y[j]-0.04),2.));
      if(spdist<sprad/3.0){
	h0_.elem(i,j) = spheight * cos(3.0*M_PI*spdist/(2.0*sprad));
	any_point_in_domain = true;
      }
    }
  }

  // Make sure that there is at least one point with nonzero initial values
  if(!any_point_in_domain){
    int i_splash = std::min(int(0.04/dx),n_boxes_-1);
    int j_splash = std::min(int(0.04/dy),n_boxes_-1);
    h0_.elem(i_splash,j_splash) = spheight;
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
  SXFunction f_step(f_step_in,f_step_out);
  f_step.init();
  cout << "generated single step dynamics (" << f_step.getAlgorithmSize() << " nodes)" << endl;
  
  // Integrate over one subinterval
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
  MXFunction f_mx(f_in,f_out);
  f_mx.init();
  cout << "generated discrete dynamics for one finite element (" << f_mx.countNodes() << " MX nodes)" << endl;
  
  // Integrate over the complete interval
  if(n_finite_elements_>1){
    f_in[0] = P;
    f_in[1] = Uk;
    f_in[2] = Vk;
    f_in[3] = Hk;
    f_inter = f_in;
    for(int j=0; j<n_finite_elements_; ++j){
      // Create a call node
      f_out = f_mx.call(f_inter);
      
      // Save intermediate state
      f_inter[1] = f_out[0];
      f_inter[2] = f_out[1];
      f_inter[3] = f_out[2];
    }
    
    // Create an integrator function
    f_mx = MXFunction(f_in,f_out);
    f_mx.init();
    cout << "generated discrete dynamics for complete interval (" << f_mx.countNodes() << " MX nodes)" << endl;    
  }
    
  // Expand the discrete dynamics
  if(false){
    SXFunction f_sx(f_mx);
    f_sx.init();
    cout << "generated discrete dynamics, SX (" << f_sx.getAlgorithmSize() << " nodes)" << endl;
    f_ = f_sx;
  } else {
    f_ = f_mx;
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


void Tester::transcribe(bool single_shooting, bool gauss_newton){
  single_shooting_ = single_shooting;
  gauss_newton_ = gauss_newton;
  
  // NLP variables
  MX P = msym("P",2);

  // Variables in the lifted NLP
  stringstream ss;
  
  // Least-squares objective function
  MX nlp_f;
  
  // Constraint function
  //  MX nlp_g = P;
  MX nlp_g = MX::sparse(0,1);
  
  // Generate full-space NLP
  MX U = u0_;  MX V = v0_;  MX H = h0_;
  for(int k=0; k<n_meas_; ++k){
    // Take a step
    MX f_arg[4] = {P,U,V,H};
    vector<MX> f_res = f_.call(vector<MX>(f_arg,f_arg+4));
    U = f_res[0];    V = f_res[1];    H = f_res[2];
    
    // Lift the height
    if(!single_shooting){
      // Initialize with measurements
      H.lift(H_meas_[k]);

      // Initialize with initial conditions
      //U.lift(u0_);
      //V.lift(v0_);
      //H.lift(DMatrix::zeros(n_boxes_  ,n_boxes_));
      
      // Initialize through simulation
      //      U.lift(U);
      //      V.lift(V);
      //      H.lift(H);
    }
    
    // Objective function term
    nlp_f.append(H-H_meas_[k]);
  }

  // Reshape to a vector
  nlp_f = flatten(nlp_f);

  // Add some regularization on p
  //  nlp_f.append(1e-3*P);

  // Form scalar objective function if not Gauss-Newton
  if(!gauss_newton_){
    nlp_f = inner_prod(nlp_f,nlp_f)/2;
  }
 
  // Function which calculates the objective terms and constraints
  vector<MX> fg_in;
  fg_in.push_back(P);

  vector<MX> fg_out;
  fg_out.push_back(nlp_f);
  fg_out.push_back(nlp_g);

  fg_mx_ = MXFunction(fg_in,fg_out);
  fg_mx_.init();
  cout << "Generated lifted NLP (" << fg_mx_.countNodes() << " nodes)" << endl;
    
}

void Tester::prepare(bool codegen, bool ipopt_as_qp_solver, bool regularization, double reg_threshold){
  verbose_ = true;
  regularization_ = regularization;
  reg_threshold_ = reg_threshold;
  
  if(codegen){
    // Make sure that command processor is available
    int flag = system(static_cast<const char*>(0));
    casadi_assert_message(flag!=0, "No command procesor available");
  }

  // Generate lifting functions
  MXFunction vdef_fcn, vinit_fcn;
  fg_mx_.generateLiftingFunctions(vdef_fcn,vinit_fcn);
  vdef_fcn.init();
  vinit_fcn.init();
  vinit_fcn_ = vinit_fcn;

  // Extract the expressions
  vector<MX> vdef_in = vdef_fcn.inputExpr();
  vector<MX> vdef_out = vdef_fcn.outputExpr();
  vector<MX> var = vdef_in;
  MX f = vdef_out.front();
  vector<MX> con(vdef_out.begin()+1,vdef_out.end());
  for(int i=0; i<con.size(); ++i){
    makeDense(con[i]);
  }

  // Get the dimensions
  nu_ = var[0].size();
  ng_ = con[0].size();
  x_.resize(var.size());
  for(int i=0; i<x_.size(); ++i){
    x_[i].n = var[i].size();
  }

  // Allocate memory
  lbu_.resize(nu_,-numeric_limits<double>::infinity());
  ubu_.resize(nu_, numeric_limits<double>::infinity());
  lbg_.resize(ng_,-numeric_limits<double>::infinity());
  ubg_.resize(ng_, numeric_limits<double>::infinity());
  g_.resize(ng_, numeric_limits<double>::quiet_NaN());
  if(!gauss_newton_){
    lambda_g_.resize(ng_,0);
    dlambda_g_.resize(ng_,0);
  }
 
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    it->init.resize(it->n,0);
    it->opt.resize(it->n,0);
    it->step.resize(it->n,0);
    if(!gauss_newton_){
      it->lam.resize(it->n,0);
      it->dlam.resize(it->n,0);
    }
  }

  // Scalar objective function
  MX obj;

  // Definition of the lifted dual variables
  vector<MX> lam_con(x_.size());

  // Multipliers
  MX lam_g;
  vector<MX> lam(x_.size());

  if(gauss_newton_){    
    // Least square objective
    obj = inner_prod(f,f)/2;
    lam_con[0] = f;
    ngL_ = lam_con[0].size();

  } else {
    // Scalar objective function
    obj = f;
    
    // Lagrange multipliers for the simple bounds on u and Lagrange multipliers corresponding to the definition of the dependent variables
    stringstream ss;
    for(int i=0; i<x_.size(); ++i){
      ss.str(string());
      ss << "lam_x" << i;
      lam[i] = msym(ss.str(),var[i].sparsity());
    }

    // Lagrange multipliers for the nonlinear constraints
    lam_g = msym("lam_g",ng_);

    if(verbose_){
      cout << "Allocated intermediate variables." << endl;
    }
   
    // Adjoint sweep to get the definitions of the lifted dual variables (Equation 3.8 in Albersmeyer2010)
    vector<vector<MX> > fseed,fsens,aseed(1),asens(1);
    aseed[0].push_back(1.0);
    aseed[0].push_back(lam_g);
    aseed[0].insert(aseed[0].end(),lam.begin()+1,lam.end());
    vdef_fcn.eval(vdef_in,vdef_out,fseed,fsens,aseed,asens,true);

    for(int i=0; i<x_.size(); ++i){
      lam_con[i] = asens[0].at(i);
    }
    ngL_ = nu_;

    if(verbose_){
      cout << "Generated the gradient of the Lagrangian." << endl;
    }
  }
  gL_.resize(ngL_, numeric_limits<double>::quiet_NaN());

  // Residual function

  // Inputs
  vector<MX> res_fcn_in;
  int n=0;
  if(!gauss_newton_){
    res_fcn_in.push_back(lam_g);       res_lam_g_ = n++;
  }
  for(int i=0; i<x_.size(); ++i){
    res_fcn_in.push_back(var[i]);      x_[i].res_var = n++;
    if(!gauss_newton_){
      res_fcn_in.push_back(lam[i]);    x_[i].res_lam = n++;
    }
  }

  // Outputs
  vector<MX> res_fcn_out;
  n=0;
  res_fcn_out.push_back(obj);               res_obj_ = n++;
  if(!gauss_newton_){
    res_fcn_out.push_back(lam_con[0]);      res_gl_ = n++;
  }
  res_fcn_out.push_back(con[0]);            res_g_ = n++;
  for(int i=1; i<x_.size(); ++i){
    res_fcn_out.push_back(con[i]-var[i]);   x_[i].res_d = n++;
    if(!gauss_newton_){
      res_fcn_out.push_back(lam_con[i]-lam[i]);  x_[i].res_lam_d = n++;
    }
  }

  // Generate function
  MXFunction res_fcn(res_fcn_in,res_fcn_out);
  res_fcn.setOption("number_of_fwd_dir",0);
  res_fcn.setOption("number_of_adj_dir",0);
  res_fcn.init();
  if(verbose_){
    cout << "Generated residual function ( " << res_fcn.getAlgorithmSize() << " nodes)." << endl;
  }

  // Generate c code and load as DLL
  if(codegen){
    dynamicCompilation(res_fcn,res_fcn_,"res_fcn","residual function");
  } else {
    res_fcn_ = res_fcn;
  }

  // Declare difference vector d and substitute out v
  vector<MX> d;
  vector<MX> d_def;
  stringstream ss;
  for(int i=1; i<x_.size(); ++i){
    ss.str(string());
    ss << "d" << i;
    MX d_i = msym(ss.str(),var[i].sparsity());
    d.push_back(d_i);
    d_def.push_back(con[i]-d_i);
  }

  // Declare difference vector lam_d and substitute out lam
  vector<MX> lam_d;
  vector<MX> lam_d_def;
  if(!gauss_newton_){
    for(int i=1; i<x_.size(); ++i){
      ss.str(string());
      ss << "lam_d" << i;
      MX lam_d_i = msym(ss.str(),var[i].sparsity());
      lam_d.push_back(lam_d_i);
      lam_d_def.push_back(lam_con[i]-lam_d_i);
    }
  }

  // Variables to be substituted and their definitions
  vector<MX> svar, sdef;
  svar.insert(svar.end(),var.begin()+1,var.end());
  sdef.insert(sdef.end(),d_def.begin(),d_def.end());
  if(!gauss_newton_){
    svar.insert(svar.end(),lam.rbegin(),lam.rend()-1);
    sdef.insert(sdef.end(),lam_d_def.rbegin(),lam_d_def.rend());    
  }

  vector<MX> ex(3);
  ex[0] = obj;
  ex[1] = con[0];
  ex[2] = lam_con[0];

  substituteInPlace(svar, sdef, ex, false);
  copy(sdef.begin(),sdef.begin()+d_def.size(),d_def.begin());  
  if(!gauss_newton_){
    copy(sdef.rbegin(),sdef.rbegin()+lam_d_def.size(),lam_d_def.begin());
  }

  MX obj_z = ex[0];
  MX g_z = ex[1];
  MX gL_z = ex[2];
  
  // Modified function Z
  vector<MX> z_in;
  n=0;
  if(!gauss_newton_){
    z_in.push_back(lam_g);                         z_lam_g_ = n++;
  }
  for(int i=0; i<x_.size(); ++i){
    z_in.push_back(i==0 ? var[0] :     d[i-1]);    x_[i].z_var = n++;
    if(!gauss_newton_){
      z_in.push_back(i==0 ? lam[0] : lam_d[i-1]);  x_[i].z_lam = n++;
    }
  }

  // Outputs
  n=0;
  vector<MX> z_out;
  z_out.push_back(obj_z);                   z_obj_ = n++;
  z_out.push_back(gL_z);                    z_gl_ = n++;
  z_out.push_back(g_z);                     z_g_ = n++;
  for(int i=1; i<x_.size(); ++i){
    z_out.push_back(d_def[i-1]);            x_[i].z_def = n++;
    if(!gauss_newton_){
      z_out.push_back(lam_d_def[i-1]);      x_[i].z_defL = n++;
    }
  }

  MXFunction zfcn(z_in,z_out);
  zfcn.setOption("name","zfcn");
  zfcn.init();
  if(verbose_){
    cout << "Generated reconstruction function ( " << zfcn.getAlgorithmSize() << " nodes)." << endl;
  }

  // Directional derivative of Z
  vector<vector<MX> > Z_fwdSeed(1,z_in);
  vector<vector<MX> > Z_fwdSens(1,z_out);
  vector<vector<MX> > Z_adjSeed;
  vector<vector<MX> > Z_adjSens;

  // Expression a + A*du in Lifted Newton (Section 2.1 in Alberspeyer2010)
  MX du = msym("du",nu_);   // Step in u
  MX dlam_g;                // Step lambda_g
  if(!gauss_newton_){
    dlam_g = msym("dlam_g",lam_g.sparsity());
  }
  
  // Interpret the Jacobian-vector multiplication as a forward directional derivative
  fill(Z_fwdSeed[0].begin(),Z_fwdSeed[0].end(),MX());
  Z_fwdSeed[0][x_[0].z_var] = du;
  for(int i=1; i<x_.size(); ++i){
    Z_fwdSeed[0][x_[i].z_var] = -d[i-1];
  }
  if(!gauss_newton_){
    Z_fwdSeed[0][z_lam_g_] = dlam_g;
    for(int i=1; i<x_.size(); ++i){
      Z_fwdSeed[0][x_[i].z_lam] = -lam_d[i-1];
    }
  }
  zfcn.eval(z_in,z_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true);    
  
  // Step expansion function inputs
  vector<MX> exp_fcn_in;
  n=0;
  if(!gauss_newton_){
    exp_fcn_in.push_back(lam_g);                         exp_lam_g_ = n++;
  }
  exp_fcn_in.push_back(du);                              exp_du_ = n++;
  if(!gauss_newton_){
    exp_fcn_in.push_back(dlam_g);                        exp_dlam_g_ = n++;
  }
  
  for(int i=0; i<x_.size(); ++i){
    exp_fcn_in.push_back(   i==0 ? var[0] :     d[i-1]);   x_[i].exp_var = n++;
    if(!gauss_newton_){
      exp_fcn_in.push_back( i==0 ? lam[0] : lam_d[i-1]);   x_[i].exp_lam = n++;
    }
  }
  
  // Step expansion function outputs
  vector<MX> exp_fcn_out;
  n=0;
  exp_fcn_out.push_back(Z_fwdSens[0][z_obj_]);           exp_osens_ = n++;
  for(int i=1; i<x_.size(); ++i){
    exp_fcn_out.push_back(Z_fwdSens[0][x_[i].z_def]);    x_[i].exp_def = n++;
    if(!gauss_newton_){
      exp_fcn_out.push_back(Z_fwdSens[0][x_[i].z_defL]); x_[i].exp_defL = n++;
    }
  }
    
  // Step expansion function
  MXFunction exp_fcn(exp_fcn_in,exp_fcn_out);
  exp_fcn.setOption("number_of_fwd_dir",0);
  exp_fcn.setOption("number_of_adj_dir",0);
  exp_fcn.setOption("name","exp_fcn");
  exp_fcn.init();
  if(verbose_){
    cout << "Generated step expansion function ( " << exp_fcn.getAlgorithmSize() << " nodes)." << endl;
  }
  
  // Generate c code and load as DLL
  if(codegen){
    dynamicCompilation(exp_fcn,exp_fcn_,"exp_fcn","step expansion function");
  } else {
    exp_fcn_ = exp_fcn;
  }
  
  // Equation (2.12 in Alberspeyer2010)
  fill(Z_fwdSeed[0].begin(),Z_fwdSeed[0].end(),MX());
  for(int i=1; i<x_.size(); ++i){
    Z_fwdSeed[0][x_[i].z_var] = d[i-1];
  }
  if(!gauss_newton_){
    for(int i=1; i<x_.size(); ++i){
      Z_fwdSeed[0][x_[i].z_lam] = lam_d[i-1];
    }
  }
  zfcn.eval(z_in,z_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true);   
  
  // Vector(s) b in Lifted Newton
  MX b_obj = z_out[z_gl_] - Z_fwdSens[0][z_gl_];
  MX b_g = z_out[z_g_] - Z_fwdSens[0][z_g_];

  // Differentiate with respect to the step to get the matrix B in Lifted Newton
  MX B_obj = zfcn.jac(x_[0].z_var,z_gl_,false,!gauss_newton_); // Exploit Hessian symmetry
  MX B_g = zfcn.jac(x_[0].z_var,z_g_);
  
  // Make sure that B_g has the right dimensions if empty
  if(B_g.empty()){
    B_g = MX::sparse(0,nu_);
  }
  
  // Get the condensed QP
  MX qpH = gauss_newton_ ? mul(trans(B_obj),B_obj) : B_obj;
  MX qpG = gauss_newton_ ? mul(trans(B_obj),b_obj) : b_obj - mul(trans(B_g),lam_g);
  MX qpA = B_g;
  MX qpB = b_g;
  
  // Make sure qpG and qpB are dense vectors
  makeDense(qpG);
  makeDense(qpB);

  if(verbose_){
    cout << "Formed qpH (" << qpH.size1() << "-by-" << qpH.size2() << ", "<< qpH.size() << " nnz)," << endl;
    cout << "       qpG (" << qpG.size1() << "-by-" << qpG.size2() << ", "<< qpG.size() << " nnz)," << endl;
    cout << "       qpA (" << qpA.size1() << "-by-" << qpA.size2() << ", "<< qpA.size() << " nnz) and" << endl;
    cout << "       qpB (" << qpB.size1() << "-by-" << qpB.size2() << ", "<< qpB.size() << " nnz)." << endl;
  }
  
  // Quadratic approximation
  vector<MX> qp_fcn_in;
  n=0;
  if(!gauss_newton_){
    qp_fcn_in.push_back(lam_g);        qpf_lam_g_ = n++;
  }
  for(int i=0; i<x_.size(); ++i){
    qp_fcn_in.push_back(var[i]);       x_[i].qpf_var = n++;
    if(!gauss_newton_){
      qp_fcn_in.push_back(lam[i]);     x_[i].qpf_lam = n++;
    }
    if(i>0){
      qp_fcn_in.push_back(d[i-1]);       x_[i].qpf_res = n++;
      if(!gauss_newton_){
	qp_fcn_in.push_back(lam_d[i-1]); x_[i].qpf_resL = n++;
      }
    }
  }
  casadi_assert(n==qp_fcn_in.size());  
  
  vector<MX> qp_fcn_out(QPF_NUM_OUT);
  qp_fcn_out[QPF_G] = qpG;
  qp_fcn_out[QPF_H] = qpH;
  qp_fcn_out[QPF_B] = qpB;
  qp_fcn_out[QPF_A] = qpA;
  
  MXFunction qp_fcn(qp_fcn_in,qp_fcn_out);
  qp_fcn.setOption("name","qp_fcn");
    qp_fcn.setOption("number_of_fwd_dir",0);
    qp_fcn.setOption("number_of_adj_dir",0);
    qp_fcn.init();
    if(verbose_){
      cout << "Generated linearization function ( " << qp_fcn.getAlgorithmSize() << " nodes)." << endl;
    }
    
    // Generate c code and load as DLL
    if(codegen){
      dynamicCompilation(qp_fcn,qp_fcn_,"qp_fcn","linearization function");
    } else {
      qp_fcn_ = qp_fcn;
    }

    // Allocate a QP solver
    if(ipopt_as_qp_solver){
      qp_solver_ = NLPQPSolver(qpH.sparsity(),qpA.sparsity());
      qp_solver_.setOption("nlp_solver",IpoptSolver::creator);
      Dictionary nlp_solver_options;
      nlp_solver_options["tol"] = 1e-12;
      nlp_solver_options["print_level"] = 0;
      nlp_solver_options["print_time"] = false;
      qp_solver_.setOption("nlp_solver_options",nlp_solver_options);
    } else {
      qp_solver_ = QPOasesSolver(qpH.sparsity(),qpA.sparsity());
      qp_solver_.setOption("printLevel","none");
    }
    
    // Initialize the QP solver
    qp_solver_.init();
    if(verbose_){
      cout << "Allocated QP solver." << endl;
    }
    
    // Residual
    for(int i=1; i<x_.size(); ++i){
      x_[i].res.resize(d[i-1].size(),0);
      if(!gauss_newton_){
	x_[i].resL.resize(lam_d[i-1].size(),0);
      }
    }
    
    if(verbose_){
      cout << "NLP preparation completed" << endl;
    }
}

void Tester::dynamicCompilation(MXFunction& f, FX& f_gen, std::string fname, std::string fdescr){
  // Append "ss" if single shooting to avoid name conflict
  if(single_shooting_){
    fname = fname + "_ss";
  }

  // C compiler
#ifdef __APPLE__
  string compiler = "gcc";
#else // __APPLE__
  string compiler = "clang";
#endif // __APPLE__

  // Optimization flag
  //string oflag = "-O3";
  string oflag = "";

  // Flag to get a DLL
#ifdef __APPLE__
  string dlflag = " -dynamiclib";
#else // __APPLE__
  string dlflag = " -shared";
#endif // __APPLE__

  // Filenames
  string cname = fname + ".c";
  string dlname = fname + ".so";
  
  // Remove existing files, if any
  string rm_command = "rm -rf " + cname + " " + dlname;
  int flag = system(rm_command.c_str());
  casadi_assert_message(flag==0, "Failed to remove old source");

  // Codegen it
  f.generateCode(cname);
  if(verbose_){
    cout << "Generated c-code for " << fdescr << " (" << cname << ")" << endl;
  }
  
  // Compile it
  string compile_command = compiler + " -fPIC " + dlflag + " " + oflag + " " + cname + " -o " + dlname;
  if(verbose_){
    cout << "Compiling " << fdescr <<  " using \"" << compile_command << "\"" << endl;
  }

  time_t time1 = time(0);
  flag = system(compile_command.c_str());
  time_t time2 = time(0);
  double comp_time = difftime(time2,time1);
  casadi_assert_message(flag==0, "Compilation failed");
  if(verbose_){
    cout << "Compiled " << fdescr << " (" << dlname << ") in " << comp_time << " s."  << endl;
  }

  // Load it
  f_gen = ExternalFunction("./" + dlname);
  f_gen.setOption("number_of_fwd_dir",0);
  f_gen.setOption("number_of_adj_dir",0);
  f_gen.setOption("name",fname + "_gen");
  f_gen.init();
  if(verbose_){
    cout << "Dynamically loaded " << fdescr << " (" << dlname << ")" << endl;
  }
}

void Tester::solve(int& iter_count){
  maxiter_ = 100;
  toldx_ = 1e-9;

  // Objective value
  obj_k_ = numeric_limits<double>::quiet_NaN();

  // Reset line-search
  meritmax_ = numeric_limits<double>::infinity();
  sigma_ = 0;

  // Current guess for the primal solution
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    copy(it->init.begin(),it->init.end(),it->opt.begin());
  }

  // Initial evaluation of the residual function
  eval_res();

  // Number of SQP iterations
  int iter = 0;

  // Reset last step-size
  pr_step_ = 0;
  du_step_ = 0;
  
  // Reset line-search
  int ls_iter = 0;
  bool ls_success = true;

  // Reset regularization
  reg_ = 0;

  // MAIN OPTIMIZATION LOOP
  while(true){
    
    // 1-norm of the primal infeasibility
    double pr_inf = primalInfeasibility();
    
    // 1-norm of the dual infeasibility
    double du_inf = dualInfeasibility();
    
    // Print header occasionally
    if(iter % 10 == 0) printIteration(cout);
    
    // Printing information about the actual iterate
    printIteration(cout,iter,obj_k_,pr_inf,du_inf,reg_,ls_iter,ls_success);

    // Checking convergence criteria
    double tol_pr_ = 1e-6; // Stopping criterion for primal infeasibility
    double tol_du_ = 1e-6; // Stopping criterion for dual infeasability
    bool converged;
    if(gauss_newton_){
      converged = iter!=0 && pr_inf < tol_pr_ && pr_step_ < toldx_; // Use step size as stopping criterion
    } else {
      converged = pr_inf < tol_pr_ && du_inf < tol_du_ && pr_step_ < toldx_; // Use lagrangian gradient as stopping criterion
    }
    if(converged){
      cout << endl;
      cout << "CasADi::SCPgen: Convergence achieved after " << iter << " iterations." << endl;
      break;
    }
    
    if (iter >= maxiter_){
      cout << endl;
      cout << "CasADi::SCPgen: Maximum number of iterations reached." << endl;
      break;
    }

    // Check if not-a-number
    if(obj_k_!=obj_k_ || pr_step_ != pr_step_ || pr_inf != pr_inf){
      cout << "CasADi::SCPgen: Aborted, nan detected" << endl;
      break;
    }
    
    // Start a new iteration
    iter++;
    
    // Form the condensed QP
    eval_qpf();

    // Regularize the QP
    if(regularization_){
      regularize();
    }

    // Solve the condensed QP
    solve_qp();

    // Expand the step
    eval_exp();
  
    // Line-search to take the step
    line_search(ls_iter, ls_success);
  }
  
  // Store optimal value
  cout << "optimal cost = " << obj_k_ << endl;
  iter_count = iter;
}

void Tester::eval_res(){
  // Pass primal variables to the residual function for initial evaluation
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    res_fcn_.setInput(it->opt,it->res_var);
  }
  
  // Pass dual variables to the residual function for initial evaluation
  if(!gauss_newton_){
    res_fcn_.setInput(lambda_g_,res_lam_g_);
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      res_fcn_.setInput(it->lam,it->res_lam);
    }
  }
  
  // Evaluate residual function
  res_fcn_.evaluate();

  // Get objective
  obj_k_ = res_fcn_.output(res_obj_).toScalar();

  // Get objective gradient
  if(!gauss_newton_){
    res_fcn_.getOutput(gL_,res_gl_);
  }

  // Get constraints
  res_fcn_.getOutput(g_,res_g_);

  // Get residuals
  for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
    res_fcn_.getOutput(it->res,  it->res_d);
    if(!gauss_newton_){
      res_fcn_.getOutput(it->resL, it->res_lam_d);
    }
  }
}

void Tester::eval_qpf(){
  // Pass primal variables to the linearization function
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    qp_fcn_.setInput(it->opt,it->qpf_var);
    if(!gauss_newton_){
      qp_fcn_.setInput(it->lam,it->qpf_lam);	
    }
  }
  
  // Pass dual variables to the linearization function
  if(!gauss_newton_){
    qp_fcn_.setInput(lambda_g_,qpf_lam_g_);
  }
  
  // Pass residual
  for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
    qp_fcn_.setInput(it->res, it->qpf_res);
    if(!gauss_newton_){
      qp_fcn_.setInput(it->resL,it->qpf_resL);
    }
  }
  
  // Evaluate to get QP
  qp_fcn_.evaluate();
}

void Tester::regularize(){
  casadi_assert(x_[0].n==2);

  DMatrix& qpH = qp_fcn_.output(QPF_H);
  
  // Regularization
  reg_ = 0;
  
  // Check the smallest eigenvalue of the Hessian
  double a = qpH.elem(0,0);
  double b = qpH.elem(0,1);
  double c = qpH.elem(1,0);
  double d = qpH.elem(1,1);
  
  // Make sure no not a numbers
  casadi_assert(a==a && b==b && c==c &&  d==d);
  
  // Make sure symmetric
  if(b!=c){
    casadi_assert_warning(fabs(b-c)<1e-10,"Hessian is not symmetric: " << b << " != " << c);
    qpH.elem(1,0) = c = b;
  }
  
  double eig_smallest = (a+d)/2 - std::sqrt(4*b*c + (a-d)*(a-d))/2;
  if(eig_smallest<reg_threshold_){
    // Regularization
    reg_ = reg_threshold_-eig_smallest;
    std::cerr << "Regularization with " << reg_ << " to ensure positive definite Hessian." << endl;
    qpH(0,0) += reg_;
    qpH(1,1) += reg_;
  }
}

void Tester::solve_qp(){
  // QP
  const DMatrix& qpH = qp_fcn_.output(QPF_H);
  const DMatrix& qpG = qp_fcn_.output(QPF_G);
  const DMatrix& qpA = qp_fcn_.output(QPF_A);
  const DMatrix& qpB = qp_fcn_.output(QPF_B);
  
  // Solve the QP
  qp_solver_.setInput(qpH,QP_H);
  qp_solver_.setInput(qpG,QP_G);
  qp_solver_.setInput(qpA,QP_A);
  std::transform(lbu_.begin(),lbu_.end(), x_[0].opt.begin(),qp_solver_.input(QP_LBX).begin(),std::minus<double>());
  std::transform(ubu_.begin(),ubu_.end(), x_[0].opt.begin(),qp_solver_.input(QP_UBX).begin(),std::minus<double>());
  
  
  //    std::transform(lbg_.begin(),lbg_.end(), qpB.begin(),qp_solver_.input(QP_LBA).begin(),std::minus<double>());
  std::transform(ubg_.begin(),ubg_.end(), qpB.begin(),qp_solver_.input(QP_UBA).begin(),std::minus<double>());
  
  //    cout << "lba = " << qp_solver_.input(QP_LBA).data() << endl;
  //    cout << "uba = " << qp_solver_.input(QP_UBA).data() << endl;
  
  qp_solver_.evaluate();
  
  // Condensed step
  const DMatrix& du = qp_solver_.output(QP_PRIMAL);
  copy(du.begin(),du.end(),x_[0].step.begin());
  
  if(!gauss_newton_){
    const DMatrix& lam_g_new = qp_solver_.output(QP_LAMBDA_A);
    copy(lam_g_new.begin(),lam_g_new.end(),dlambda_g_.begin());
    std::transform(dlambda_g_.begin(),dlambda_g_.end(),lambda_g_.begin(),dlambda_g_.begin(),std::minus<double>());
    
    const DMatrix& dlam_u = qp_solver_.output(QP_LAMBDA_X);
    copy(dlam_u.begin(),dlam_u.end(),x_[0].dlam.begin());
  }  
}

void Tester::line_search(int& ls_iter, bool& ls_success){
  // Calculate penalty parameter of merit function
  sigma_ = std::max(sigma_,1.01*norm_inf(qp_solver_.output(QP_LAMBDA_X).data()));
  sigma_ = std::max(sigma_,1.01*norm_inf(qp_solver_.output(QP_LAMBDA_A).data()));
  
  // Calculate L1-merit function in the actual iterate
  double l1_infeas = primalInfeasibility();
  
  // Right-hand side of Armijo condition
  double F_sens = exp_fcn_.output(exp_osens_).toScalar();
  double L1dir = F_sens - sigma_ * l1_infeas;
  double L1merit = obj_k_ + sigma_ * l1_infeas;
  
  // Store the actual merit function
  if(L1merit<meritmax_){
    meritmax_ = L1merit;
  }
  
  // Stepsize
  double t = 1.0, t_prev = 0.0;
  double fk_cand;
  
  // Merit function value in candidate
  double L1merit_cand = 0;
  
  // Reset line-search counter, success marker
  ls_iter = 0;
  ls_success = false;
  int maxiter_ls_ = 1;
  
  // Line-search parameter, restoration factor of stepsize
  double beta_ = 0.8;
  
  // Armijo condition, coefficient of decrease in merit
  double c1_ = 1e-4;
  
  // Line-search
  //log("Starting line-search");
  
  // Line-search loop
  while (true){
    
    // Take the primal step
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      for(int i=0; i<it->n; ++i){
	it->opt[i] += (t-t_prev) * it->step[i];
      }
    }
    
    // Take the dual step
    if(!gauss_newton_){
      for(int i=0; i<ng_; ++i){
	lambda_g_[i] += (t-t_prev) * dlambda_g_[i];
      }
      for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
	if(it==x_.begin()){
	  //copy(it->dlam.begin(),it->dlam.end(),it->lam.begin()); // BUG?
	} else {
	  for(int i=0; i<it->n; ++i){
	    it->lam[i] += (t-t_prev) * it->dlam[i];
	  }
	}
      }
    }
    
    // Evaluate residual function to get objective and constraints (and residuals for the next iteration)
    eval_res();
    ls_iter++;      
    
    // Calculating merit-function in candidate
    l1_infeas = primalInfeasibility();
    L1merit_cand = obj_k_ + sigma_ * l1_infeas;
    
    // Calculating maximal merit function value so far
    if (L1merit_cand <= meritmax_ + t * c1_ * L1dir){
      
      // Accepting candidate
      ls_success = true;
      //log("Line-search completed, candidate accepted");
      break;
    }
    
    // Line-search not successful, but we accept it.
    if(ls_iter == maxiter_ls_){
      //log("Line-search completed, maximum number of iterations");
      break;
    }
    
    // Backtracking
    t_prev = t;
    t = beta_ * t;
  }

  // Calculate primal step-size
  pr_step_ = 0;
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    for(vector<double>::const_iterator i=it->step.begin(); i!=it->step.end(); ++i){
      pr_step_ += fabs(*i);
    }
  }
  pr_step_ *= t;

  // Calculate the dual step-size
  if(!gauss_newton_){
    du_step_ = 0;
    for(vector<double>::const_iterator i=dlambda_g_.begin(); i!=dlambda_g_.end(); ++i){
      du_step_ += fabs(*i);
    }

    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      if(it==x_.begin()){
	// ?
      } else {
	for(vector<double>::const_iterator i=it->dlam.begin(); i!=it->dlam.end(); ++i){
	  du_step_ += fabs(*i);
	}
      }
    }
    du_step_ *= t;
  }
}

void Tester::eval_exp(){
  // Pass primal step/variables
  exp_fcn_.setInput(x_[0].step, exp_du_);
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    if(it==x_.begin()){
      exp_fcn_.setInput(it->opt,it->exp_var);
    } else {
      exp_fcn_.setInput(it->res,it->exp_var);
    }
  }
  
  // Pass dual step/variables
  if(!gauss_newton_){
    exp_fcn_.setInput(dlambda_g_,exp_dlam_g_);
    exp_fcn_.setInput(lambda_g_,exp_lam_g_);
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      if(it==x_.begin()){
	exp_fcn_.setInput(it->lam,it->exp_lam);
      } else {
	exp_fcn_.setInput(it->resL,it->exp_lam);
      }
    }
  }

  // Perform the step expansion
  exp_fcn_.evaluate();

  // Expanded primal step
  for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
    const DMatrix& dv = exp_fcn_.output(it->exp_def);
    copy(dv.begin(),dv.end(),it->step.begin());
  }
  
  // Expanded dual step
  if(!gauss_newton_){
    for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
      const DMatrix& dlam_v = exp_fcn_.output(it->exp_defL);
      copy(dlam_v.begin(),dlam_v.end(),it->dlam.begin());
    }
  }  
}

void Tester::optimize(double drag_guess, double depth_guess, int& iter_count, double& sol_time, double& drag_est, double& depth_est){
  if(verbose_){
    cout << "Starting parameter estimation" << endl;
  }
  
  // Initial guess
  x_[0].init[0] = drag_guess;
  x_[0].init[1] = depth_guess;

  if(x_.size()>1){
    // Initialize lifted variables using the generated function
    vinit_fcn_.setInput(x_[0].init);
    vinit_fcn_.evaluate();    
    for(int i=1; i<x_.size(); ++i){
      vinit_fcn_.getOutput(x_[i].init,i-1);
    }
  }
  if(verbose_){
    cout << "Passed initial guess" << endl;
  }

  // Reset dual guess
  if(!gauss_newton_){
    fill(lambda_g_.begin(),lambda_g_.end(),0);
    fill(dlambda_g_.begin(),dlambda_g_.end(),0);
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      fill(it->lam.begin(),it->lam.end(),0);
      fill(it->dlam.begin(),it->dlam.end(),0);
    }
  }

  // Bounds on the variables
  lbu_.at(0) = 1e-4; // drag positive
  lbu_.at(1) = 1e-4; // depth positive

  ubu_.at(0) = 1e2; // max drag
  ubu_.at(1) = 1e2; // max depth

  // Constraint bounds
  //  fill(lbg_.begin(),lbg_.end(),0);

    fill(ubg_.begin(),ubg_.end(),10);


  clock_t time1 = clock();
  solve(iter_count);
  clock_t time2 = clock();
  
  // Solution statistics  
  sol_time = double(time2-time1)/CLOCKS_PER_SEC;
  drag_est = x_[0].opt.at(0);
  depth_est = x_[0].opt.at(1);
}

double Tester::primalInfeasibility(){
  // L1-norm of the primal infeasibility
  double pr_inf = 0;
  
  // Variable bounds
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    if(it==x_.begin()){
      
      // Simple bounds
      for(int i=0; i<it->n; ++i) pr_inf +=  ::fmax(it->opt[i]-ubu_[i],0.);
      for(int i=0; i<it->n; ++i) pr_inf +=  ::fmax(lbu_[i]-it->opt[i],0.);

    } else {
      
      // Lifted variables
      for(int i=0; i<it->n; ++i) pr_inf += ::fabs(it->res[i]);
    }
  }
  
  // Nonlinear bounds
  for(int i=0; i<ng_; ++i) pr_inf += ::fmax(g_[i]-ubg_[i],0.);
  for(int i=0; i<ng_; ++i) pr_inf += ::fmax(lbg_[i]-g_[i],0.);
  
  return pr_inf;
}  

double Tester::dualInfeasibility(){
  // Not implemented for Gauss-Newton
  if(gauss_newton_) return 0;

  // L1-norm of the dual infeasibility
  double du_inf = 0;
  
  // Lifted variables
  for(int i=0; i<ngL_; ++i) du_inf += ::fabs(gL_[i]);

  return du_inf;
}

void Tester::printIteration(std::ostream &stream){
  stream << setw(4)  << "iter";
  stream << setw(14) << "objective";
  stream << setw(11) << "inf_pr";
  stream << setw(11) << "inf_du";
  stream << setw(11) << "pr_step";
  stream << setw(11) << "du_step";
  stream << setw(8) << "lg(rg)";
  stream << setw(3) << "ls";
  stream << ' ';

  // Problem specific: generic functionality needed
  stream << setw(9) << "drag";
  stream << setw(9) << "depth";

  stream << endl;
}
  
void Tester::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double rg, int ls_trials, bool ls_success){
  stream << setw(4) << iter;
  stream << scientific;
  stream << setw(14) << setprecision(6) << obj;
  stream << setw(11) << setprecision(2) << pr_inf;
  stream << setw(11);
  if(gauss_newton_){
    stream << "-";
  } else {
    stream << setprecision(2) << du_inf;
  }
  stream << setw(11) << setprecision(2) << pr_step_;
  stream << setw(11);
  if(gauss_newton_){
    stream << "-";
  } else {
    stream << setprecision(2) << du_step_;
  }
  stream << fixed;
  if(rg>0){
    stream << setw(8) << setprecision(2) << log10(rg);
  } else {
    stream << setw(8) << "-";
  }
  stream << setw(3) << ls_trials;
  stream << (ls_success ? ' ' : 'F');

  // Problem specific: generic functionality needed
  stream << setw(9) << setprecision(3) << x_[0].opt.at(0);
  stream << setw(9) << setprecision(3) << x_[0].opt.at(1);

  stream << endl;
}

int main(){

  // True parameter values
  double drag_true = 2.0, depth_true = 0.01;
  
  // Use IPOPT as QP solver (can handle non-convex QPs)
  bool ipopt_as_qp_solver = true;

  // Use Gauss-Newton method
  bool gauss_newton = false;

  // Codegen the Lifted Newton functions
  bool codegen = false;
  
  // Regularize the QP
  bool regularization = false;

  // Smallest allowed eigenvalue for the regularization
  double reg_threshold = 1e-8;

  // Problem size
  // int  n_boxes = 100, n_euler = 100, n_finite_elements = 1, n_meas = 20;
  //int  n_boxes = 30, n_euler = 40, n_finite_elements = 25, n_meas = 20; // Paper
  int n_boxes = 15, n_euler = 20, n_finite_elements = 1, n_meas = 20;

  // Initial guesses
  vector<double> drag_guess, depth_guess;
  drag_guess.push_back( 2.0); depth_guess.push_back(0.01); // Optimal solution
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
  Tester t(n_boxes,n_euler,n_finite_elements,n_meas);
    
  // Perform the modelling
  t.model();

  // Optimization parameters
  t.simulate(drag_true, depth_true);
  
  // For both single and multiple shooting
  for(int sol=0; sol<2; ++sol){

    // Transcribe as an NLP
    bool single_shooting = sol==0;
    t.transcribe(single_shooting,gauss_newton);
  
    // Prepare the NLP solver
    t.prepare(codegen, ipopt_as_qp_solver, regularization, reg_threshold);

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

