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
  enum LinOut{LIN_QPG,LIN_H,LIN_QPB,LIN_QPA,LIN_NUM_OUT};
  
  /// Generate initial guess for lifted variables
  FX vinit_fcn_;

  /// Residual function
  FX rfcn_;
  
  /// Quadratic approximation
  FX lfcn_;
  
  /// Step expansion
  FX efcn_;
  
  /// Dimensions
  int nu_, ng_, ngL_;

  // Simple and nonlinear bounds
  vector<double> lbu_, ubu_, lbg_, ubg_;

  /// Multipliers for the nonlinear bounds
  vector<double> lambda_g_, dlambda_g_;

  int g_con_;
  int g_f_, g_g_;

  int z_con_;
  int z_obj_;
  int z_g_;

  int l_con_;

  int e_du_, e_dlam_g_, e_con_;

  struct Var{
    int n;

    int g_var, g_lam;
    int g_d, g_lam_d;

    int z_var, z_lam;
    int z_def, z_defL;

    int l_var, l_lam, l_res, l_resL;

    int e_var, e_lam;
    int e_exp, e_expL;

    vector<double> step, init, opt, lam, dlam;
    vector<double> res, resL;
    
  };
  vector<Var> x_;  
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
  vector<MX> var = vdef_fcn.inputExpr();
  vector<MX> con = vdef_fcn.outputExpr();

  MX f = con.front();
  con.erase(con.begin());
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

  // Lagrangian gradient of the objective function with respect to u and v
  vector<MX> gL(x_.size());

  // Multipliers
  MX lam_g;
  vector<MX> lam(x_.size());

  if(gauss_newton_){
    
    // Least square objective
    obj = inner_prod(f,f)/2;
    gL[0] = f;
    ngL_ = gL[0].size();

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
   
    // Lagrangian
    MX lag = obj;
    if(!con[0].empty()) lag += inner_prod(lam_g,con[0]);
    if(!var[0].empty()) lag += inner_prod(lam[0],var[0]);
    for(int i=1; i<x_.size(); ++i){
      if(!con[i].empty()) lag += inner_prod(flatten(lam[i]),flatten(con[i])); // note: equal to trace(mul(lam.T,con))
    }

    // Lagrangian function
    MXFunction lfcn(var,lag);
    lfcn.init();

    // Adjoint sweep to get gradient
    vector<MX> lfcn_in = lfcn.inputExpr();
    vector<MX> lfcn_out = lfcn.outputExpr();
    vector<vector<MX> > fseed,fsens,aseed(1),asens(1);
    aseed[0].resize(1,1.0);
    asens[0].resize(2);
    lfcn.eval(lfcn_in,lfcn_out,fseed,fsens,aseed,asens,true);
    for(int i=0; i<x_.size(); ++i){
      gL[i] = asens[0].at(i);
      makeDense(gL[i]);
    }
    ngL_ = nu_;

    if(verbose_){
      cout << "Generated the gradient of the Lagrangian." << endl;
    }
  }

  // Residual function

  // Inputs
  vector<MX> rfcn_in;
  int n=0;
  if(!gauss_newton_){
    rfcn_in.push_back(lam_g);       g_con_ = n++;
  }
  for(int i=0; i<x_.size(); ++i){
    rfcn_in.push_back(var[i]);      x_[i].g_var = n++;
    if(!gauss_newton_){
      rfcn_in.push_back(lam[i]);    x_[i].g_lam = n++;
    }
  }

  // Outputs
  vector<MX> rfcn_out;
  n=0;
  rfcn_out.push_back(obj);               g_f_ = n++;
  rfcn_out.push_back(con[0]);            g_g_ = n++;
  for(int i=1; i<x_.size(); ++i){
    rfcn_out.push_back(con[i]-var[i]);   x_[i].g_d = n++;
    if(!gauss_newton_){
      rfcn_out.push_back(gL[i]-lam[i]);  x_[i].g_lam_d = n++;
    }
  }

  // Generate function
  MXFunction rfcn(rfcn_in,rfcn_out);
  rfcn.setOption("number_of_fwd_dir",0);
  rfcn.setOption("number_of_adj_dir",0);
  rfcn.init();
  if(verbose_){
    cout << "Generated residual function ( " << rfcn.getAlgorithmSize() << " nodes)." << endl;
  }

  // Generate c code and load as DLL
  if(codegen){
    dynamicCompilation(rfcn,rfcn_,"rfcn","residual function");
  } else {
    rfcn_ = rfcn;
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

  // Declare difference vector lam_d and substitute out lam_v
  vector<MX> lam_d;
  vector<MX> lam_d_def;
  if(!gauss_newton_){
    for(int i=1; i<x_.size(); ++i){
      ss.str(string());
      ss << "lam_d" << i;
      MX lam_d_i = msym(ss.str(),var[i].sparsity());
      lam_d.push_back(lam_d_i);
      lam_d_def.push_back(lam[i]-lam_d_i);
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

  vector<MX> ex(2);
  ex[0] = obj;
  ex[1] = con[0];
  ex.insert(ex.end(),gL.begin(),gL.end());
 
  substituteInPlace(svar, sdef, ex, false);
  copy(sdef.begin(),sdef.begin()+d_def.size(),d_def.begin());  
  if(!gauss_newton_){
    copy(sdef.rbegin(),sdef.rbegin()+lam_d_def.size(),lam_d_def.begin());
  }

  MX obj_z = ex[0];
  MX g_z = ex[1];
  vector<MX> gL_z(ex.begin()+2,ex.end());
  
  // Modified function Z
  vector<MX> z_in;
  n=0;
  if(!gauss_newton_){
    z_in.push_back(lam_g);                       z_con_ = n++;
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
  z_out.push_back(gL_z[0]);                 z_obj_ = n++;
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

  // Step expansion function
  if(x_.size()>1){
    vector<MX> v_exp(x_.size());
    vector<MX> vL_exp(x_.size());
    
    // Step in u and lambda_g
    MX du = msym("du",nu_);
    MX dlam_g;
    if(!gauss_newton_){
      dlam_g = msym("dlam_g",lam_g.sparsity());
    }

    Z_fwdSeed[0][x_[0].z_var] = du;
    for(int i=1; i<x_.size(); ++i){
      Z_fwdSeed[0][x_[i].z_var] = -d[i-1];
    }

    if(!gauss_newton_){
      Z_fwdSeed[0][z_con_] = dlam_g;
      Z_fwdSeed[0][x_[0].z_lam] = MX();
      for(int i=1; i<x_.size(); ++i){
	Z_fwdSeed[0][x_[i].z_lam] = -lam_d[i-1];
      }
    }

    zfcn.eval(z_in,z_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true);
    
    for(int i=1; i<x_.size(); ++i){
      v_exp[i] = Z_fwdSens[0][x_[i].z_def];
      if(!gauss_newton_){
	vL_exp[i] = Z_fwdSens[0][x_[i].z_defL];
      }
    }
    
    // Function input
    vector<MX> efcn_in;
    n=0;
    if(!gauss_newton_){
      efcn_in.push_back(lam_g);                         e_con_ = n++;
    }
    efcn_in.push_back(du);                              e_du_ = n++;
    if(!gauss_newton_){
      efcn_in.push_back(dlam_g);                        e_dlam_g_ = n++;
    }
    
    for(int i=0; i<x_.size(); ++i){
      efcn_in.push_back(   i==0 ? var[0] :     d[i-1]);   x_[i].e_var = n++;
      if(!gauss_newton_){
	efcn_in.push_back( i==0 ? lam[0] : lam_d[i-1]);   x_[i].e_lam = n++;
      }
    }

    // Function outputs
    vector<MX> efcn_out;
    n=0;
    for(int i=1; i<x_.size(); ++i){
      efcn_out.push_back(v_exp[i]);                     x_[i].e_exp = n++;
      if(!gauss_newton_){
	efcn_out.push_back(vL_exp[i]);                  x_[i].e_expL = n++;
      }
    }
    
    // Form function object
    MXFunction efcn(efcn_in,efcn_out);
    efcn.setOption("number_of_fwd_dir",0);
    efcn.setOption("number_of_adj_dir",0);
    efcn.setOption("name","efcn");
    efcn.init();
    if(verbose_){
      cout << "Generated step expansion function ( " << efcn.getAlgorithmSize() << " nodes)." << endl;
    }

    // Generate c code and load as DLL
    if(codegen){
      dynamicCompilation(efcn,efcn_,"efcn","step expansion function");
    } else {
      efcn_ = efcn;
    }
  }
  
  // Matrix A and B in lifted Newton  
  MX qpH = zfcn.jac(x_[0].z_var,z_obj_,false,!gauss_newton_); // Exploit Hessian symmetry
  MX qpA = zfcn.jac(x_[0].z_var,z_g_,false);

  // Make sure that A has the right dimensions if empty
  if(qpA.empty()){
    qpA = MX::sparse(0,nu_);
  }

  if(verbose_){
    cout << "Formed qpH (dimension " << qpH.size1() << "-by-" << qpH.size2() << ", "<< qpH.size() << " nonzeros) " <<
               "and qpA (dimension " << qpA.size1() << "-by-" << qpA.size2() << ", "<< qpA.size() << " nonzeros)." << endl;
  }
  
  MX qpG = gL_z[0];
  MX qpB = g_z;

  if(x_.size()>1){
    Z_fwdSeed[0][x_[0].z_var] = MX();
    for(int i=1; i<x_.size(); ++i){
      Z_fwdSeed[0][x_[i].z_var] = -d[i-1];
    }
    
    if(!gauss_newton_){
      Z_fwdSeed[0][z_con_] = MX();
      Z_fwdSeed[0][x_[0].z_lam] = MX();
      for(int i=1; i<x_.size(); ++i){
	Z_fwdSeed[0][x_[i].z_lam] = MX();
      }
    }
    
    zfcn.eval(z_in,z_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true);
    
    qpG += Z_fwdSens[0][z_obj_];
    qpB += Z_fwdSens[0][z_g_];
  }
  if(verbose_){
    cout << "Formed qpG (dimension " << qpG.size1() << "-by-" << qpG.size2() << ", "<< qpG.size() << " nonzeros) " <<
               "and qpB (dimension " << qpB.size1() << "-by-" << qpB.size2() << ", "<< qpB.size() << " nonzeros)." << endl;
  }
  
  // Generate Gauss-Newton Hessian
  if(gauss_newton_){
    qpG = mul(trans(qpH),qpG);
    qpH = mul(trans(qpH),qpH);
    if(verbose_){
      cout << "Gauss Newton Hessian (dimension " << qpH.size1() << "-by-" << qpH.size2() << ", "<< qpH.size() << " nonzeros)." << endl;
    }
  }
  
  // Make sure qpG and qpB are dense vectors
  makeDense(qpG);
  makeDense(qpB);

  // Quadratic approximation
  vector<MX> lfcn_in;
  n=0;
  if(!gauss_newton_){
    lfcn_in.push_back(lam_g);        l_con_ = n++;
  }
  for(int i=0; i<x_.size(); ++i){
    lfcn_in.push_back(var[i]);       x_[i].l_var = n++;
    if(!gauss_newton_){
      lfcn_in.push_back(lam[i]);     x_[i].l_lam = n++;
    }
    if(i>0){
      lfcn_in.push_back(d[i-1]);       x_[i].l_res = n++;
      if(!gauss_newton_){
	lfcn_in.push_back(lam_d[i-1]); x_[i].l_resL = n++;
      }
    }
  }
  casadi_assert(n==lfcn_in.size());  
  
  vector<MX> lfcn_out(LIN_NUM_OUT);
  lfcn_out[LIN_QPG] = qpG;
  lfcn_out[LIN_H] = qpH;
  lfcn_out[LIN_QPB] = qpB;
  lfcn_out[LIN_QPA] = qpA;
  

  MXFunction lfcn(lfcn_in,lfcn_out);
  lfcn.setOption("name","lfcn");
  lfcn.setOption("number_of_fwd_dir",0);
  lfcn.setOption("number_of_adj_dir",0);
  lfcn.init();
  if(verbose_){
    cout << "Generated linearization function ( " << lfcn.getAlgorithmSize() << " nodes)." << endl;
  }

  // Generate c code and load as DLL
  if(codegen){
    dynamicCompilation(lfcn,lfcn_,"lfcn","linearization function");
  } else {
    lfcn_ = lfcn;
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
  double f_k = numeric_limits<double>::quiet_NaN();
  
  // Current guess for the primal solution
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    copy(it->init.begin(),it->init.end(),it->opt.begin());
  }
  
  int k=0;  
  while(true){

    // Pass primal variables to the residual function
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      rfcn_.setInput(it->opt,it->g_var);
    }

    // Pass dual variables to the residual function
    if(!gauss_newton_){
      rfcn_.setInput(lambda_g_,g_con_);
      for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
	rfcn_.setInput(it->lam,it->g_lam);
      }
    }

    // Evaluate residual function
    rfcn_.evaluate();

    // Store residual
    f_k = rfcn_.output(g_f_).toScalar();
    const DMatrix& r_g = rfcn_.output(g_g_);
    for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
      rfcn_.getOutput(it->res,  it->g_d);
      if(!gauss_newton_){
	rfcn_.getOutput(it->resL, it->g_lam_d);
      }
    }
    
    // Pass primal variables to the linearization function
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      lfcn_.setInput(it->opt,it->l_var);
      if(!gauss_newton_){
	lfcn_.setInput(it->lam,it->l_lam);	
      }
    }
    
    // Pass dual variables to the linearization function
    if(!gauss_newton_){
      lfcn_.setInput(lambda_g_,l_con_);
    }

    // Pass residual
    for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
      lfcn_.setInput(it->res, it->l_res);
      if(!gauss_newton_){
	lfcn_.setInput(it->resL,it->l_resL);
      }
    }

    // Evaluate to get QP
    lfcn_.evaluate();
    
    // QP
    DMatrix& qpH = lfcn_.output(LIN_H);
    const DMatrix& qpG = lfcn_.output(LIN_QPG);
    const DMatrix& qpA = lfcn_.output(LIN_QPA);
    const DMatrix& qpB = lfcn_.output(LIN_QPB);

    // Regularization
    double reg = 0;
    
    // Check the smallest eigenvalue of the Hessian
    if(regularization_ && x_[0].n==2){
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
	reg = reg_threshold_-eig_smallest;
	std::cerr << "Regularization with " << reg << " to ensure positive definite Hessian." << endl;
	qpH(0,0) += reg;
	qpH(1,1) += reg;
      }
    }
    
    
    // Solve the QP
    qp_solver_.setInput(qpH,QP_H);
    qp_solver_.setInput(qpG,QP_G);
    qp_solver_.setInput(qpA,QP_A);
    std::transform(lbu_.begin(),lbu_.end(), x_[0].opt.begin(),qp_solver_.input(QP_LBX).begin(),std::minus<double>());
    std::transform(ubu_.begin(),ubu_.end(), x_[0].opt.begin(),qp_solver_.input(QP_UBX).begin(),std::minus<double>());
    std::transform(lbg_.begin(),lbg_.end(), qpB.begin(),qp_solver_.input(QP_LBA).begin(),std::minus<double>());
    std::transform(ubg_.begin(),ubg_.end(), qpB.begin(),qp_solver_.input(QP_UBA).begin(),std::minus<double>());
    qp_solver_.evaluate();
    const DMatrix& du = qp_solver_.output(QP_PRIMAL);
    const DMatrix& dlam_u = qp_solver_.output(QP_LAMBDA_X);
    const DMatrix& dlam_g = qp_solver_.output(QP_LAMBDA_A);    
    
    // Condensed step
    copy(du.begin(),du.end(),x_[0].step.begin());
    if(!gauss_newton_){
      copy(dlam_g.begin(),dlam_g.end(),dlambda_g_.begin());
      copy(dlam_u.begin(),dlam_u.end(),x_[0].dlam.begin());
    }

    // Expand the step
    if(x_.size()>1){
      efcn_.setInput(du,e_du_);
      efcn_.setInput(lfcn_.input(x_[0].l_var),x_[0].e_var);
      for(int i=1; i<x_.size(); ++i){
	efcn_.setInput(lfcn_.input(x_[i].l_res),x_[i].e_var);
      }
      if(!gauss_newton_){
	efcn_.setInput(dlam_g,e_dlam_g_);
	efcn_.setInput(lfcn_.input(l_con_),e_con_);
	efcn_.setInput(lfcn_.input(x_[0].l_lam),x_[0].e_lam);
	for(int i=1; i<x_.size(); ++i){
	  efcn_.setInput(lfcn_.input(x_[i].l_resL),x_[i].e_lam);
	}
      }
      efcn_.evaluate();
      
      // Expanded primal step
      for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
	const DMatrix& dv = efcn_.output(it->e_exp);
	copy(dv.begin(),dv.end(),it->step.begin());
      }
      
      // Expanded dual step
      if(!gauss_newton_){
	for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
	  const DMatrix& dlam_v = efcn_.output(it->e_expL);
	  copy(dlam_v.begin(),dlam_v.end(),it->dlam.begin());
	}
      }
    }
    
    // Take a full step
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      transform(it->step.begin(),it->step.end(),it->opt.begin(),it->opt.begin(),plus<double>());
    }

    if(!gauss_newton_){
      transform(dlambda_g_.begin(),dlambda_g_.end(),lambda_g_.begin(),lambda_g_.begin(),plus<double>());
      for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
	//	if(it==x_.begin()){
	  copy(it->dlam.begin(),it->dlam.end(),it->lam.begin()); // BUG?
	  //	} else {
	  transform(it->dlam.begin(),it->dlam.end(),it->lam.begin(),it->lam.begin(),plus<double>()); // BUG?
	  //	}
      }
    }

    // Step size
    double norm_step=0;
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      for(vector<double>::const_iterator i=it->step.begin(); i!=it->step.end(); ++i){
	norm_step += *i**i;
      }
    }
    if(!gauss_newton_){
      for(vector<double>::const_iterator i=dlambda_g_.begin();	i!=dlambda_g_.end(); ++i){
	norm_step += *i**i;
      }
      for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
	for(vector<double>::const_iterator i=it->dlam.begin(); i!=it->dlam.end(); ++i){
	  norm_step += *i**i;
	}
      }
    }
    norm_step = sqrt(norm_step);
    
    // Constraint violation
    double norm_viol = 0;

    // Variable bounds
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      if(it==x_.begin()){

	// Simple bounds
	for(int i=0; i<it->n; ++i){
	  double d = ::fmax(it->opt[i]-ubu_[i],0.) + ::fmax(lbu_[i]-it->opt[i],0.);
	  norm_viol += d*d;
	}
      } else {

	// Lifted variables
	for(int i=0; i<it->n; ++i){
	  double d = it->res[i];
	  norm_viol += d*d;
	}
      }
    }

    // Nonlinear bounds
    for(int i=0; i<ng_; ++i){
      double d = ::fmax(r_g.at(i)-ubg_[i],0.) + ::fmax(lbg_[i]-r_g.at(i),0.);
      norm_viol += d*d;
    }
        
    // Square root to get 2-norm
    norm_viol = sqrt(norm_viol);
    
    // Print progress (including the header every 10 rows)
    if(k % 10 == 0){
      cout << setw(4) << "iter" << setw(20) << "objective" << setw(20) << "drag"          << setw(20) << "depth"          << setw(20) << "norm_step" << setw(20) << "norm_viol" << endl;
    }
    cout   << setw(4) <<     k  << setw(20) <<  f_k        << setw(20) << x_[0].opt.at(0) << setw(20) << x_[0].opt.at(1) << setw(20) <<  norm_step  << setw(20) <<  norm_viol  << endl;
    
    
    // Check if stopping criteria is satisfied
    if(norm_viol + norm_step < toldx_){
      cout << "Convergence achieved!" << endl;
      break;
    }

    // Check if not-a-number
    if(f_k!=f_k || norm_step != norm_step || norm_viol != norm_viol){
      cout << "Aborted, nan detected" << endl;
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
  fill(lbg_.begin(),lbg_.end(),0);
  fill(ubg_.begin(),ubg_.end(),0);

  clock_t time1 = clock();
  solve(iter_count);
  clock_t time2 = clock();
  
  // Solution statistics  
  sol_time = double(time2-time1)/CLOCKS_PER_SEC;
  drag_est = x_[0].opt.at(0);
  depth_est = x_[0].opt.at(1);
}

int main(){

  // True parameter values
  double drag_true = 2.0, depth_true = 0.01;
  
  // Use IPOPT as QP solver (can handle non-convex QPs)
  bool ipopt_as_qp_solver = false;

  // Use Gauss-Newton method
  bool gauss_newton = true;

  // Codegen the Lifted Newton functions
  bool codegen = false;
  
  // Regularize the QP
  bool regularization = true;

  // Smallest allowed eigenvalue for the regularization
  double reg_threshold = 1e-8;

  // Problem size
  // int  n_boxes = 100, n_euler = 100, n_finite_elements = 1, n_meas = 20;
  // int  n_boxes = 30, n_euler = 40, n_finite_elements = 25, n_meas = 20; // Paper
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

