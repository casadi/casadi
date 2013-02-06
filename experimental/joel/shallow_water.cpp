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
  void transcribe(bool single_shooting, bool gauss_newton);

  // Prepare NLP
  void prepare();

  // Solve NLP
  void solve(int& iter_count);

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

  // Options
  bool single_shooting_;
  bool verbose_;
  bool gauss_newton_;

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

  /// Residual
  vector<double> d_, lam_d_;
   
  /// Primal step
  vector<double> du_,dv_;
   
  /// Dual step
  vector<double> dlambda_u_, dlambda_h_, dlambda_g_;
  
  /// Indices
  enum GIn{G_U,G_LAM_U,G_LAM_G,G_V,G_LAM_H,G_NUM_IN};
  enum GOut{G_D,G_LAM_D,G_G,G_F,G_NUM_OUT};
  
  enum LinIn{LIN_U,LIN_LAM_U,LIN_LAM_G,LIN_V,LIN_LAM_H,LIN_D,LIN_LAM_D,LIN_NUM_IN};
  enum LinOut{LIN_F1,LIN_J1,LIN_G,LIN_J2,LIN_NUM_OUT};
  
  enum ExpIn{EXP_U,EXP_LAM_U,EXP_LAM_G,EXP_D,EXP_LAM_D,EXP_DU,EXP_DLAM_G,EXP_NUM_IN};
  enum ExpOut{EXP_V,EXP_LAM_H,EXP_NUM_OUT};

  enum ZIn{Z_U,Z_D,Z_LAM_D,Z_LAM_U,Z_LAM_G,Z_NUM_IN};
  enum ZOut{Z_D_DEF,Z_LAM_D_DEF,Z_FG,Z_NUM_OUT};

  /// Residual function
  FX rfcn_;
  
  /// Quadratic approximation
  FX lfcn_;
  
  /// Step expansion
  FX efcn_;
  
  /// Dimensions
  int nu_, nv_, ng_, ngL_;
  
  vector<double> u_init_, lbu_, ubu_, u_opt_, lambda_u_;
  vector<double> lbg_, ubg_, lambda_g_;
  vector<double> v_init_, lbv_, ubv_, v_opt_;
  vector<double> lambda_h_;
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
  MXFunction f_mx(f_in,f_out);
  f_mx.init();
  cout << "generated discrete dynamics, MX (" << f_mx.countNodes() << " nodes)" << endl;
  
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
  MX nlp_u = msym("u",2);
  MX nlp_v = msym("v",single_shooting ? 0 : n_boxes_*n_boxes_*n_meas_);

  // Variables in the lifted NLP
  MX P = nlp_u;
  int v_offset = 0;
  
  // Least-squares objective function
  MX nlp_f;
  
  // Constraint function
  MX nlp_g = MX::sparse(0,1);
  
  // Lifted equations
  MX nlp_h = MX::sparse(0,1);

  // Generate full-space NLP
  MX U = u0_;  MX V = v0_;  MX H = h0_;
  for(int k=0; k<n_meas_; ++k){
    // Take a step
    MX f_arg[4] = {P,U,V,H};
    vector<MX> f_res = f_.call(vector<MX>(f_arg,f_arg+4));
    U = f_res[0];    V = f_res[1];    H = f_res[2];
    
    if(!single_shooting){
      // Lift the variable
      MX H_def = H;
      H = nlp_v[Slice(v_offset,v_offset+n_boxes_*n_boxes_)];
      H = reshape(H,h0_.sparsity());
      v_offset += H.size();
      
      // Constraint function term
      nlp_h.append(flatten(H_def));
    }
    
    // Objective function term
    nlp_f.append(flatten(H-H_meas_[k]));
  }

  // Form scalar objective function if not Gauss-Newton
  if(!gauss_newton_){
    nlp_f = inner_prod(nlp_f,nlp_f)/2;
  }
 
  // Function which calculates the objective terms and constraints
  vector<MX> fg_in;
  fg_in.push_back(nlp_u);
  fg_in.push_back(nlp_v);

  vector<MX> fg_out;
  fg_out.push_back(nlp_f);
  fg_out.push_back(nlp_g);
  fg_out.push_back(nlp_h);

  fg_mx_ = MXFunction(fg_in,fg_out);
  fg_mx_.init();
  cout << "Generated lifted NLP (" << fg_mx_.countNodes() << " nodes)" << endl;
  
  // Expand NLP
  fg_sx_ = SXFunction(fg_mx_);
  fg_sx_.init();
  cout << "expanded lifted NLP (" << fg_sx_.getAlgorithmSize() << " nodes)" << endl;
  
  // Prepare the NLP solver
  prepare();
}

void Tester::prepare(){
  verbose_ = true;

  // Generate lifting functions
  //  MXFunction F,G,Z;
  //  fg_mx_.generateLiftingFunctions(F,G,Z);
  //  F.init();
  //  G.init();
  //  Z.init();

  // Extract the free variables and split into independent and dependent variables
  SXMatrix u = fg_sx_.inputExpr(0);
  SXMatrix v = fg_sx_.inputExpr(1);

  // Extract the constraint equations and split into constraints and definitions of dependent variables
  SXMatrix f = fg_sx_.outputExpr(0);
  SXMatrix g = fg_sx_.outputExpr(1);
  SXMatrix h = fg_sx_.outputExpr(2);
  makeDense(h);

  // Get the dimensions
  nu_ = u.size();
  nv_ = v.size();
  ng_ = g.size();
 
  u_init_.resize(nu_,0);
  u_opt_.resize(nu_,0);
  lbu_.resize(nu_,-numeric_limits<double>::infinity());
  ubu_.resize(nu_, numeric_limits<double>::infinity());
  lambda_u_.resize(nu_,0);

  v_init_.resize(nv_,0);
  v_opt_.resize(nv_,0);
  lbv_.resize(nv_,-numeric_limits<double>::infinity());
  ubv_.resize(nv_, numeric_limits<double>::infinity());

  lbg_.resize(ng_,-numeric_limits<double>::infinity());
  ubg_.resize(ng_, numeric_limits<double>::infinity());
  lambda_g_.resize(ng_,0);
  lambda_h_.resize(nv_,0);
        
  // Scalar objective function
  SXMatrix obj;
  
  // Lagrangian gradient of the objective function with respect to u and v
  SXMatrix gL_u, gL_v;

  // Multipliers
  SXMatrix lam_u, lam_g, lam_h;
  if(gauss_newton_){
    
    // Least square objective
    obj = inner_prod(f,f)/2;
    gL_u = f;
    ngL_ = gL_u.size();

  } else {
    
    // Scalar objective function
    obj = f;
    
    // Lagrange multipliers for the simple bounds on u
    lam_u = ssym("lam_u",nu_);
        
    // Lagrange multipliers corresponding to the definition of the dependent variables
    lam_h = ssym("lam_h",nv_);

    // Lagrange multipliers for the nonlinear constraints that aren't eliminated
    lam_g = ssym("lam_g",ng_);

    if(verbose_){
      cout << "Allocated intermediate variables." << endl;
    }
    
    // Lagrangian
    SXMatrix lag = obj;
    if(!u.empty()) lag += inner_prod(lam_u,u);
    if(!g.empty()) lag += inner_prod(lam_g,g);
    if(!h.empty()) lag += inner_prod(lam_h,h);

    // Lagrangian function
    vector<SXMatrix> x(2);
    x[0] = u;
    x[1] = v;
    SXFunction lfcn(x,lag);
    lfcn.init();

    // Adjoint sweep to get gradient
    vector<SXMatrix> lfcn_in = lfcn.inputExpr();
    vector<SXMatrix> lfcn_out = lfcn.outputExpr();
    vector<vector<SXMatrix> > fseed,fsens,aseed(1),asens(1);
    aseed[0].resize(1,1.0);
    asens[0].resize(2);
    lfcn.eval(lfcn_in,lfcn_out,fseed,fsens,aseed,asens,true);
    gL_u = asens[0].at(0);
    gL_v = asens[0].at(1);
    makeDense(gL_u);
    makeDense(gL_v);
    ngL_ = nu_;

    if(verbose_){
      cout << "Generated the gradient of the Lagrangian." << endl;
    }
  }

  // Residual function G
  vector<SXMatrix> G_in(G_NUM_IN);
  G_in[G_U] = u;
  G_in[G_V] = v;
  G_in[G_LAM_U] = lam_u;
  G_in[G_LAM_G] = lam_g;
  G_in[G_LAM_H] = lam_h;

  vector<SXMatrix> G_out(G_NUM_OUT);
  G_out[G_D] = h-v;
  G_out[G_LAM_D] = gL_v-lam_h;
  G_out[G_G] = vertcat(h-v,g);
  G_out[G_F] = obj;
  rfcn_ = SXFunction(G_in,G_out);
  rfcn_.setOption("number_of_fwd_dir",0);
  rfcn_.setOption("number_of_adj_dir",0);
  rfcn_.init();
  if(verbose_){
    cout << "Generated residual function ( " << shared_cast<SXFunction>(rfcn_).getAlgorithmSize() << " nodes)." << endl;
  }

  // Declare difference vector d and substitute out v
  SXMatrix d = ssym("d",nv_);
  SXMatrix d_def = h-d;
  vector<SXMatrix> ex(4);
  ex[0] = gL_u;
  ex[1] = g;
  ex[2] = obj;
  ex[3] = gL_v;
  substituteInPlace(v, d_def, ex, false);
  SXMatrix gL_u_z = ex[0];
  SXMatrix g_z = ex[1];
  SXMatrix obj_z = ex[2];
  SXMatrix gL_v_z = ex[3];

  // Declare difference vector lam_d and substitute out lam_h
  SXMatrix lam_d;
  SXMatrix lam_d_def;
  if(!gauss_newton_){
    lam_d = ssym("lam_d",nv_);
    lam_d_def = lam_h-lam_d;
    ex.resize(2);
    substituteInPlace(lam_h, lam_d_def, ex, false);    
    gL_u_z = ex[0];
    g_z = ex[1];
  }
  
  // Modified function Z
  vector<SXMatrix> zfcn_in(Z_NUM_IN);
  zfcn_in[Z_U] = u;
  zfcn_in[Z_D] = d;
  zfcn_in[Z_LAM_D] = lam_d;
  zfcn_in[Z_LAM_U] = lam_u;
  zfcn_in[Z_LAM_G] = lam_g;
  
  vector<SXMatrix> zfcn_out(Z_NUM_OUT);
  zfcn_out[Z_D_DEF] = d_def;
  zfcn_out[Z_LAM_D_DEF] = lam_d_def;
  zfcn_out[Z_FG] = vertcat(gL_u_z,g_z);
  
  SXFunction zfcn(zfcn_in,zfcn_out);
  zfcn.setOption("name","zfcn");
  zfcn.init();
  if(verbose_){
    cout << "Generated reconstruction function ( " << zfcn.getAlgorithmSize() << " nodes)." << endl;
  }

  // Matrix A and B in lifted Newton
  SXMatrix B = zfcn.jac(Z_U,Z_FG);
  SXMatrix B1 = B(Slice(0,ngL_),Slice(0,B.size2()));
  SXMatrix B2 = B(Slice(ngL_,B.size1()),Slice(0,B.size2()));
  if(verbose_){
    cout << "Formed B1 (dimension " << B1.size1() << "-by-" << B1.size2() << ", "<< B1.size() << " nonzeros) " <<
    "and B2 (dimension " << B2.size1() << "-by-" << B2.size2() << ", "<< B2.size() << " nonzeros)." << endl;
  }
  
  // Step in u
  SXMatrix du = ssym("du",nu_);
  SXMatrix dlam_g = ssym("dlam_g",lam_g.sparsity());
  
  SXMatrix b1 = gL_u_z;
  SXMatrix b2 = g_z;
  SXMatrix e, eL;
  if(nv_ > 0){
    
    // Directional derivative of Z
    vector<vector<SXMatrix> > Z_fwdSeed(2,zfcn_in);
    vector<vector<SXMatrix> > Z_fwdSens(2,zfcn_out);
    vector<vector<SXMatrix> > Z_adjSeed;
    vector<vector<SXMatrix> > Z_adjSens;
    
    Z_fwdSeed[0][Z_U].setZero();
    Z_fwdSeed[0][Z_D] = -d;
    Z_fwdSeed[0][Z_LAM_U].setZero();
    Z_fwdSeed[0][Z_LAM_G].setZero();
    
    Z_fwdSeed[1][Z_U] = du;
    Z_fwdSeed[1][Z_D] = -d;
    Z_fwdSeed[1][Z_LAM_U].setZero();
    Z_fwdSeed[1][Z_LAM_G] = dlam_g;
    
    zfcn.eval(zfcn_in,zfcn_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true);
    
    b1 += Z_fwdSens[0][Z_FG](Slice(0,ngL_));
    b2 += Z_fwdSens[0][Z_FG](Slice(ngL_,B.size1()));
    e = Z_fwdSens[1][Z_D_DEF];
    eL = Z_fwdSens[1][Z_LAM_D_DEF];
  }
  if(verbose_){
    cout << "Formed b1 (dimension " << b1.size1() << "-by-" << b1.size2() << ", "<< b1.size() << " nonzeros) " <<
    "and b2 (dimension " << b2.size1() << "-by-" << b2.size2() << ", "<< b2.size() << " nonzeros)." << endl;
  }
  
  // Generate Gauss-Newton Hessian
  if(gauss_newton_){
    b1 = mul(trans(B1),b1);
    B1 = mul(trans(B1),B1);
    if(verbose_){
      cout << "Gauss Newton Hessian (dimension " << B1.size1() << "-by-" << B1.size2() << ", "<< B1.size() << " nonzeros)." << endl;
    }
  }
  
  // Make sure b1 and b2 are dense vectors
  makeDense(b1);
  makeDense(b2);
  
  // Quadratic approximation
  vector<SXMatrix> lfcn_in(LIN_NUM_IN);
  lfcn_in[LIN_U] = u;
  lfcn_in[LIN_V] = v;
  lfcn_in[LIN_D] = d;
  lfcn_in[LIN_LAM_D] = lam_d;
  lfcn_in[LIN_LAM_U] = lam_u;
  lfcn_in[LIN_LAM_G] = lam_g;
  lfcn_in[LIN_LAM_H] = lam_h;
  
  vector<SXMatrix> lfcn_out(LIN_NUM_OUT);
  lfcn_out[LIN_F1] = b1;
  lfcn_out[LIN_J1] = B1;
  lfcn_out[LIN_G] = b2;
  lfcn_out[LIN_J2] = B2;
  lfcn_ = SXFunction(lfcn_in,lfcn_out);
  lfcn_.setOption("name","lfcn");
  lfcn_.setOption("number_of_fwd_dir",0);
  lfcn_.setOption("number_of_adj_dir",0);
  lfcn_.init();
  if(verbose_){
    cout << "Generated linearization function ( " << shared_cast<SXFunction>(lfcn_).getAlgorithmSize() << " nodes)." << endl;
  }
    
  // Step expansion
  vector<SXMatrix> efcn_in(EXP_NUM_IN);
  efcn_in[EXP_U] = u;
  efcn_in[EXP_D] = d;
  efcn_in[EXP_LAM_D] = lam_d;
  efcn_in[EXP_LAM_U] = lam_u;
  efcn_in[EXP_LAM_G] = lam_g;
  efcn_in[EXP_DU] = du;
  efcn_in[EXP_DLAM_G] = dlam_g;
  vector<SXMatrix> efcn_out(EXP_NUM_OUT);
  efcn_out[EXP_V] = e;
  efcn_out[EXP_LAM_H] = eL;
  efcn_ = SXFunction(efcn_in,efcn_out);
  efcn_.setOption("number_of_fwd_dir",0);
  efcn_.setOption("number_of_adj_dir",0);
  efcn_.setOption("name","efcn");
  efcn_.init();
  if(verbose_){
    cout << "Generated step expansion function ( " << shared_cast<SXFunction>(efcn_).getAlgorithmSize() << " nodes)." << endl;
  }
  
  // Allocate a QP solver
  qp_solver_ = QPOasesSolver(B1.sparsity(),B2.sparsity());
  qp_solver_.setOption("printLevel","none");
  
  // Initialize the QP solver
  qp_solver_.init();
  if(verbose_){
    cout << "Allocated QP solver." << endl;
  }

  // Residual
  d_.resize(d.size(),0);
  lam_d_.resize(lam_d.size(),0);
  
  // Primal step
  du_.resize(nu_);
  dv_.resize(nv_);

  // Dual step
  dlambda_u_.resize(nu_);
  dlambda_g_.resize(ng_);
  dlambda_h_.resize(nv_);
}

void Tester::solve(int& iter_count){
  maxiter_ = 100;
  toldx_ = 1e-9;

  // Objective value
  double f_k = numeric_limits<double>::quiet_NaN();
  
  // Current guess for the primal solution
  copy(u_init_.begin(),u_init_.end(),u_opt_.begin());
  copy(v_init_.begin(),v_init_.end(),v_opt_.begin());
  
  int k=0;
  
  // Does G depend on the multipliers?
  bool has_lam_u =  !gauss_newton_ && nu_>0;
  bool has_lam_g =  !gauss_newton_ && ng_>0;
  bool has_lam_h =  !gauss_newton_ && nv_>0;
  
  while(true){

    // Evaluate residual
    rfcn_.setInput(u_opt_,G_U);
    rfcn_.setInput(v_opt_,G_V);

    if(has_lam_u) rfcn_.setInput(lambda_u_,G_LAM_U);
    if(has_lam_g) rfcn_.setInput(lambda_g_,G_LAM_G);
    if(has_lam_h) rfcn_.setInput(lambda_h_,G_LAM_H);

    rfcn_.evaluate();
    rfcn_.getOutput(d_,G_D);
    rfcn_.getOutput(lam_d_,G_LAM_D);
    f_k = rfcn_.output(G_F).toScalar();
    const DMatrix& g_k = rfcn_.output(G_G);
    
    // Construct the QP
    lfcn_.setInput(u_opt_,LIN_U);
    lfcn_.setInput(v_opt_,LIN_V);
    
    if(has_lam_u) lfcn_.setInput(lambda_u_,LIN_LAM_U);
    if(has_lam_g) lfcn_.setInput(lambda_g_,LIN_LAM_G);
    if(has_lam_h) lfcn_.setInput(lambda_h_,LIN_LAM_H);

    lfcn_.setInput(d_,LIN_D);
    lfcn_.setInput(lam_d_,LIN_LAM_D);
    lfcn_.evaluate();
    DMatrix& B1_k = lfcn_.output(LIN_J1);
    const DMatrix& b1_k = lfcn_.output(LIN_F1);
    const DMatrix& B2_k = lfcn_.output(LIN_J2);
    const DMatrix& b2_k = lfcn_.output(LIN_G);

    // Regularization
    double reg = 0;
    bool regularization = true;
    
    // Check the smallest eigenvalue of the Hessian
    if(regularization && nu_==2){
      double a = B1_k.elem(0,0);
      double b = B1_k.elem(0,1);
      double c = B1_k.elem(1,0);
      double d = B1_k.elem(1,1);
      
      // Make sure no not a numbers
      casadi_assert(a==a && b==b && c==c &&  d==d);
      
      // Make sure symmetric
      if(b!=c){
	casadi_assert_warning(fabs(b-c)<1e-10,"Hessian is not symmetric: " << b << " != " << c);
	B1_k.elem(1,0) = c = b;
      }
      
      double eig_smallest = (a+d)/2 - std::sqrt(4*b*c + (a-d)*(a-d))/2;
      double threshold = 1e-8;
      if(eig_smallest<threshold){
	// Regularization
	reg = threshold-eig_smallest;
	std::cerr << "Regularization with " << reg << " to ensure positive definite Hessian." << endl;
	B1_k(0,0) += reg;
	B1_k(1,1) += reg;
      }
    }
    
    
    // Solve the QP
    qp_solver_.setInput(B1_k,QP_H);
    qp_solver_.setInput(b1_k,QP_G);
    qp_solver_.setInput(B2_k,QP_A);
    std::transform(lbu_.begin(),lbu_.end(),u_opt_.begin(),qp_solver_.input(QP_LBX).begin(),std::minus<double>());
    std::transform(ubu_.begin(),ubu_.end(),u_opt_.begin(),qp_solver_.input(QP_UBX).begin(),std::minus<double>());
    std::transform(lbg_.begin(),lbg_.end(), b2_k.begin(),qp_solver_.input(QP_LBA).begin(),std::minus<double>());
    std::transform(ubg_.begin(),ubg_.end(), b2_k.begin(),qp_solver_.input(QP_UBA).begin(),std::minus<double>());
    qp_solver_.evaluate();
    const DMatrix& du = qp_solver_.output(QP_PRIMAL);
    const DMatrix& dlam_u = qp_solver_.output(QP_LAMBDA_X);
    const DMatrix& dlam_g = qp_solver_.output(QP_LAMBDA_A);    
    
    // Expand the step
    efcn_.setInput(lfcn_.input(LIN_U),EXP_U);
    efcn_.setInput(lfcn_.input(LIN_LAM_U),EXP_LAM_U);
    efcn_.setInput(lfcn_.input(LIN_LAM_G),EXP_LAM_G);
    efcn_.setInput(lfcn_.input(LIN_D),EXP_D);
    efcn_.setInput(lfcn_.input(LIN_LAM_D),EXP_LAM_D);
    efcn_.setInput(du,EXP_DU);
    if(has_lam_g) efcn_.setInput(dlam_g,EXP_DLAM_G);
    efcn_.evaluate();
    const DMatrix& dv = efcn_.output(EXP_V);
    const DMatrix& dlam_h = efcn_.output(EXP_LAM_H);
    
    // Expanded primal step
    copy(du.begin(),du.end(),du_.begin());
    copy(dv.begin(),dv.end(),dv_.begin());

    // Expanded dual step
    copy(dlam_u.begin(),dlam_u.end(),dlambda_u_.begin());
    copy(dlam_g.begin(),dlam_g.end(),dlambda_g_.begin());
    copy(dlam_h.begin(),dlam_h.end(),dlambda_h_.begin());
    
    // Take a full step
    transform(du_.begin(),du_.end(),u_opt_.begin(),u_opt_.begin(),plus<double>());
    transform(dv_.begin(),dv_.end(),v_opt_.begin(),v_opt_.begin(),plus<double>());

    copy(dlambda_u_.begin(),dlambda_u_.end(),lambda_u_.begin());
    transform(dlambda_g_.begin(),dlambda_g_.end(),lambda_g_.begin(),lambda_g_.begin(),plus<double>());
    transform(dlambda_h_.begin(),dlambda_h_.end(),lambda_h_.begin(),lambda_h_.begin(),plus<double>());

    // Step size
    double norm_step=0;
    for(vector<double>::const_iterator it=du_.begin(); it!=du_.end(); ++it)  norm_step += *it**it;
    for(vector<double>::const_iterator it=dv_.begin(); it!=dv_.end(); ++it)  norm_step += *it**it;
    if(!gauss_newton_){
      for(vector<double>::const_iterator it=dlambda_g_.begin(); it!=dlambda_g_.end(); ++it) norm_step += *it**it;
      for(vector<double>::const_iterator it=dlambda_h_.begin(); it!=dlambda_h_.end(); ++it) norm_step += *it**it;
    }
    norm_step = sqrt(norm_step);
    
    // Constraint violation
    double norm_viol = 0;
    for(int i=0; i<nu_; ++i){
      double d = ::fmax(u_opt_.at(i)-ubu_.at(i),0.) + ::fmax(lbu_.at(i)-u_opt_.at(i),0.);
      norm_viol += d*d;
    }
    for(int i=0; i<nv_; ++i){
      double d = ::fmax(v_opt_.at(i)-ubv_.at(i),0.) + ::fmax(lbv_.at(i)-v_opt_.at(i),0.);
      norm_viol += d*d;
    }
    for(int i=0; i<ng_; ++i){
      double d = ::fmax(g_k.at(nv_+i)-ubg_.at(i),0.) + ::fmax(lbg_.at(i)-g_k.at(nv_+i),0.);
      norm_viol += d*d;
    }
    for(int i=0; i<nv_; ++i){
      double d = ::fmax(g_k.at(i),0.) + ::fmax(-g_k.at(i),0.);
      norm_viol += d*d;
    }

    norm_viol = sqrt(norm_viol);
    
    // Print progress (including the header every 10 rows)
    if(k % 10 == 0){
      cout << setw(4) << "iter" << setw(20) << "objective" << setw(20) << "norm_step" << setw(20) << "norm_viol" << endl;
    }
    cout   << setw(4) <<     k  << setw(20) <<  f_k        << setw(20) <<  norm_step  << setw(20) <<  norm_viol  << endl;
    
    
    // Check if stopping criteria is satisfied
    if(norm_viol + norm_step < toldx_){
      cout << "Convergence achieved!" << endl;
      break;
    }
    
    // Increase iteration count
    k = k+1;
    
    // Check if number of iterations have been reached
    if(k >= maxiter_){
      cout << "Maximum number of iterations (" << maxiter_ << ") reached" << endl;
      break;
    }
  }
  
  // Store optimal value
  cout << "optimal cost = " << f_k << endl;
  iter_count = k;
}

void Tester::optimize(double drag_guess, double depth_guess, int& iter_count, double& sol_time, double& drag_est, double& depth_est){
  // Initial guess for the parameters
  u_init_[0] = drag_guess;
  u_init_[1] = depth_guess;
  fill(v_init_.begin(),v_init_.end(),0);

  // Initial guess for the heights
  if(!single_shooting_){
    vector<double>::iterator it=v_init_.begin();
    for(int k=0; k<n_meas_; ++k){
      copy(H_meas_[k].begin(),H_meas_[k].end(),it);
      it += n_boxes_*n_boxes_;
    }
  }
  
  fill(lbg_.begin(),lbg_.end(),0);
  fill(ubg_.begin(),ubg_.end(),0);
  lbu_[0] = 0;
  lbu_[1] = 0;

  clock_t time1 = clock();
  solve(iter_count);
  clock_t time2 = clock();
  
  // Solution statistics  
  sol_time = double(time2-time1)/CLOCKS_PER_SEC;
  drag_est = u_opt_.at(0);
  depth_est = u_opt_.at(1);
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
    bool gauss_newton = true;
    t.transcribe(single_shooting,gauss_newton);
  
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

