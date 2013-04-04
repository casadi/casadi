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

#include "scpgen_internal.hpp"
#include "symbolic/casadi.hpp"
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>

#ifdef WITH_DL 
#include <cstdlib>
#endif // WITH_DL 

using namespace std;
namespace CasADi{

SCPgenInternal::SCPgenInternal(const FX& F, const FX& G, const FX& H, const FX& J) : NLPSolverInternal(F,G,H,J){
  casadi_warning("SCPgen is under development");
  addOption("qp_solver",         OT_QPSOLVER,   GenericType(),    "The QP solver to be used by the SQP method");
  addOption("qp_solver_options", OT_DICTIONARY, GenericType(),    "Options to be passed to the QP solver");
  addOption("hessian_approximation", OT_STRING, "limited-memory", "limited-memory|exact");
  addOption("maxiter",           OT_INTEGER,      50,             "Maximum number of SQP iterations");
  addOption("maxiter_ls",        OT_INTEGER,       1,             "Maximum number of linesearch iterations");
  addOption("tol_pr",            OT_REAL,       1e-6,             "Stopping criterion for primal infeasibility");
  addOption("tol_du",            OT_REAL,       1e-6,             "Stopping criterion for dual infeasability");
  addOption("tol_reg",           OT_REAL,       1e-11,            "Stopping criterion for regularization");
  addOption("tol_pr_step",       OT_REAL,       1e-6,             "Stopping criterion for the step size");
  addOption("c1",                OT_REAL,       1e-4,             "Armijo condition, coefficient of decrease in merit");
  addOption("beta",              OT_REAL,       0.8,              "Line-search parameter, restoration factor of stepsize");
  addOption("merit_memsize",     OT_INTEGER,      4,              "Size of memory to store history of merit function values");
  addOption("merit_start",       OT_REAL,      1e-8,              "Lower bound for the merit function parameter");
  addOption("lbfgs_memory",      OT_INTEGER,     10,              "Size of L-BFGS memory.");
  addOption("regularize",        OT_BOOLEAN,  false,              "Automatic regularization of Lagrange Hessian.");
  addOption("print_header",      OT_BOOLEAN,   true,              "Print the header with problem statistics");
  addOption("codegen",           OT_BOOLEAN,  false,              "C-code generation");
  addOption("reg_threshold",     OT_REAL,      1e-8,              "Threshold for the regularization.");
  addOption("name_x",      OT_STRINGVECTOR,  GenericType(),       "Names of the variables.");
  addOption("print_x",           OT_INTEGERVECTOR,  GenericType(), "Which variables to print.");
  addOption("compiler",          OT_STRING,    "gcc -fPIC -O2",    "Compiler command to be used for compiling generated code");
  addOption("print_time",        OT_BOOLEAN, true,                 "Print information about execution time");
  
  // Monitors
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx", true);
}


SCPgenInternal::~SCPgenInternal(){
}

void SCPgenInternal::init(){
  // Call the init method of the base class
  NLPSolverInternal::init();
    
  // Read options
  maxiter_ = getOption("maxiter");
  maxiter_ls_ = getOption("maxiter_ls");
  c1_ = getOption("c1");
  beta_ = getOption("beta");
  lbfgs_memory_ = getOption("lbfgs_memory");
  tol_pr_ = getOption("tol_pr");
  tol_du_ = getOption("tol_du");
  tol_reg_ = getOption("tol_reg");
  regularize_ = getOption("regularize");
  codegen_ = getOption("codegen");
  reg_threshold_ = getOption("reg_threshold");
  print_time_ = getOption("print_time");
  tol_pr_step_ = getOption("tol_pr_step");
  merit_memsize_ = getOption("merit_memsize");
  merit_start_ = getOption("merit_start");

  // Name the components
  if(hasSetOption("name_x")){
    name_x_ = getOption("name_x");
    casadi_assert(name_x_.size()==nx_);
  } else {
    stringstream ss;
    name_x_.resize(nx_);
    for(int i=0; i<nx_; ++i){
      ss.str(string());
      ss << "x" << i;
      name_x_[i] = ss.str();
    }
  }

  // Components to print
  if(hasSetOption("print_x")){
    print_x_ = getOption("print_x");
  } else {
    print_x_.resize(0);
  }

  //  if(getOption("hessian_approximation")=="exact")
  //    hess_mode_ = HESS_EXACT;
  //  else if(getOption("hessian_approximation")=="limited-memory")
  //    hess_mode_ = HESS_BFGS;
   
  //  if (hess_mode_== HESS_EXACT && H_.isNull()) {
  //    if (!getOption("generate_hessian")){
  //      casadi_error("SCPgenInternal::evaluate: you set option 'hessian_approximation' to 'exact', but no hessian was supplied. Try with option \"generate_hessian\".");
  //   }
  //  }
  
  //  // If the Hessian is generated, we use exact approximation by default
  //  if (bool(getOption("generate_hessian"))){
  //    setOption("hessian_approximation", "exact");
  //  }

  // Form a function that calculates noth the objective and constraints
  casadi_assert(!F_.isNull() || !G_.isNull());
  casadi_assert(F_.isNull() || is_a<MXFunction>(F_));
  casadi_assert(G_.isNull() || is_a<MXFunction>(G_));
  MXFunction F = shared_cast<MXFunction>(F_);
  MXFunction G = shared_cast<MXFunction>(G_);
  
  // Get the expressions to be able to assemble the fg function below
  MX nlp_x, nlp_p = msym("p",0,1), nlp_f, nlp_g;
  if(F.isNull()){
    // Root-finding problem
    nlp_x = G.inputExpr(0);
    nlp_g = G.outputExpr(0);
    if(G.getNumInputs()>1){
      nlp_p = G.inputExpr(1);
    }
  } else {
    // Minimization problem
    nlp_x = F.inputExpr(0);
    nlp_f = F.outputExpr(0);
    if(F.getNumInputs()>1){
      nlp_p = F.inputExpr(1);
    }
    if(G.isNull()){
      nlp_g = MX::sparse(0,1);
    } else {
      nlp_g = G.outputExpr(0);
      if(G.getNumInputs()>F.getNumInputs()){
	nlp_p = G.inputExpr(1);
      }
    }
  }
  
  vector<MX> fg_in;
  fg_in.push_back(nlp_x);
  fg_in.push_back(nlp_p);

  vector<MX> fg_out;
  fg_out.push_back(nlp_f);
  fg_out.push_back(nlp_g);

  MXFunction fg(fg_in,fg_out);
  fg.init();
  
  if(codegen_){
#ifdef WITH_DL 
    // Make sure that command processor is available
    int flag = system(static_cast<const char*>(0));
    casadi_assert_message(flag!=0, "No command procesor available");
#else // WITH_DL 
    casadi_error("Codegen in SCPgen requires CasADi to be compiled with option \"WITH_DL\" enabled");
#endif // WITH_DL 
  }

  // Generate lifting functions
  MXFunction vdef_fcn, vinit_fcn;
  fg.generateLiftingFunctions(vdef_fcn,vinit_fcn);
  vdef_fcn.init();
  vinit_fcn.init();
  vinit_fcn_ = vinit_fcn;

  // Extract the expressions
  vector<MX> vdef_in = vdef_fcn.inputExpr();
  vector<MX> vdef_out = vdef_fcn.outputExpr();

  // Get the dimensions
  x_ = vdef_in.at(0);
  p_ = vdef_in.at(1);
  v_.resize(vdef_in.size()-2);
  for(int i=0; i<v_.size(); ++i){
    v_[i].v = vdef_in.at(i+2);
    v_[i].v_def = vdef_out.at(i+2);
    v_[i].n = v_[i].v.size();
  }

  // Allocate memory
  x_lb_.resize(nx_,-numeric_limits<double>::infinity());
  x_ub_.resize(nx_, numeric_limits<double>::infinity());
  g_lb_.resize(ng_,-numeric_limits<double>::infinity());
  g_ub_.resize(ng_, numeric_limits<double>::infinity());
  g_.resize(ng_, numeric_limits<double>::quiet_NaN());
  if(!gauss_newton_){
    g_lam_.resize(ng_,0);
    g_dlam_.resize(ng_,0);
  }
  qpH_times_du_.resize(nx_);
 
  x_init_.resize(nx_,0);
  x_opt_.resize(nx_,0);
  x_step_.resize(nx_,0);
  if(!gauss_newton_){
    x_lam_.resize(nx_,0);
    x_dlam_.resize(nx_,0);
  }

  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    it->init.resize(it->n,0);
    it->opt.resize(it->n,0);
    it->step.resize(it->n,0);
    if(!gauss_newton_){
      it->lam.resize(it->n,0);
      it->dlam.resize(it->n,0);
    }
  }

  // Line-search memory
  merit_mem_.resize(merit_memsize_);

  // Scalar objective function
  MX obj;

  // Multipliers
  MX g_lam, x_lam;

  // Definition of the lifted dual variables
  MX p_defL, gL_defL;

  if(gauss_newton_){    
    // Least square objective
    obj = inner_prod(vdef_out[0],vdef_out[0])/2;
    gL_defL = vdef_out[0];
    ngL_ = gL_defL.size();

  } else {
    // Scalar objective function
    obj = vdef_out[0];
    
    // Lagrange multipliers for the simple bounds on u and Lagrange multipliers corresponding to the definition of the dependent variables
    x_lam = msym("x_lam",x_.sparsity());
    stringstream ss;
    int i=0;
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      ss.str(string());
      ss << "lam_x" << i++;
      it->v_lam = msym(ss.str(),it->v.sparsity());
    }
    
    // Lagrange multipliers for the nonlinear constraints
    g_lam = msym("g_lam",ng_);

    if(verbose_){
      cout << "Allocated intermediate variables." << endl;
    }
   
    // Adjoint sweep to get the definitions of the lifted dual variables (Equation 3.8 in Albersmeyer2010)
    vector<vector<MX> > fseed,fsens,aseed(1),asens(1);
    aseed[0].push_back(1.0);
    aseed[0].push_back(g_lam);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      aseed[0].push_back(it->v_lam);
    }
    vdef_fcn.eval(vdef_in,vdef_out,fseed,fsens,aseed,asens,true);
    i=0;

    gL_defL = asens[0].at(i++);
    if(gL_defL.isNull()) gL_defL = MX(x_.sparsity()); // Needed?

    p_defL = asens[0].at(i++);
    if(p_defL.isNull()) p_defL = MX(p_.sparsity()); // Needed?

    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      it->v_defL = asens[0].at(i++);
      if(it->v_defL.isNull()){
        it->v_defL = MX(it->v.sparsity());
      }
    }

    // Multipliers for the simple bounds
    gL_defL += x_lam;
    ngL_ = nx_;

    if(verbose_){
      cout << "Generated the gradient of the Lagrangian." << endl;
    }
  }
  gL_.resize(ngL_, numeric_limits<double>::quiet_NaN());

  // Residual function

  // Inputs
  vector<MX> res_fcn_in;
  int n=0;
  res_fcn_in.push_back(x_);             res_x_ = n++;
  res_fcn_in.push_back(p_);             res_p_ = n++;
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    res_fcn_in.push_back(it->v);        it->res_var = n++;
  }
  if(!gauss_newton_){
    res_fcn_in.push_back(g_lam);        res_g_lam_ = n++;
    res_fcn_in.push_back(x_lam);        res_x_lam_ = n++;
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      res_fcn_in.push_back(it->v_lam);  it->res_lam = n++;
    }
  }

  // Outputs
  vector<MX> res_fcn_out;
  n=0;
  res_fcn_out.push_back(obj);                            res_obj_ = n++;
  res_fcn_out.push_back(gL_defL);                        res_gl_ = n++;
  res_fcn_out.push_back(vdef_out[1]);                    res_g_ = n++;
  res_fcn_out.push_back(p_defL);                         res_p_d_ = n++;
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    res_fcn_out.push_back(it->v_def - it->v);             it->res_d = n++;
  }
    
  if(!gauss_newton_){
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      res_fcn_out.push_back(it->v_defL - it->v_lam);     it->res_lam_d = n++;
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
  if(codegen_){
    dynamicCompilation(res_fcn,res_fcn_,"res_fcn","residual function");
  } else {
    res_fcn_ = res_fcn;
  }

  // Declare difference vector d and substitute out p and v
  stringstream ss;
  int i=0;
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    ss.str(string());
    ss << "d" << i++;
    it->d = msym(ss.str(),it->v.sparsity());
    it->d_def = it->v_def - it->d;
  }

  // Declare difference vector lam_d and substitute out lam
  if(!gauss_newton_){
    int i=0;
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      ss.str(string());
      ss << "d_lam" << i++;
      it->d_lam = msym(ss.str(),it->v.sparsity());
      it->d_defL = it->v_defL - it->d_lam;
    }
  }

  // Variables to be substituted and their definitions
  vector<MX> svar, sdef;
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    svar.push_back(it->v);
    sdef.push_back(it->d_def);
  }
  if(!gauss_newton_){
    for(vector<Var>::reverse_iterator it=v_.rbegin(); it!=v_.rend(); ++it){
      svar.push_back(it->v_lam);
      sdef.push_back(it->d_defL);
    }
  }

  vector<MX> ex(4);
  ex[0] = obj;
  ex[1] = vdef_out[1];
  ex[2] = gL_defL;
  ex[3] = p_defL;

  substituteInPlace(svar, sdef, ex, false);
  i=0;
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    it->d_def = sdef[i++];
  }
  if(!gauss_newton_){
    for(vector<Var>::reverse_iterator it=v_.rbegin(); it!=v_.rend(); ++it){
      it->d_defL = sdef[i++];
    }
  }

  MX obj_z = ex[0];
  MX g_z = ex[1];
  MX gL_z = ex[2];
  MX p_z = ex[3];
  
  // Modified function inputs
  vector<MX> mfcn_in;
  n=0;
  mfcn_in.push_back(p_);                               mod_p_ = n++;
  mfcn_in.push_back(x_);                               mod_x_ = n++;
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    mfcn_in.push_back(it->d);                          it->mod_var = n++;
  }

  // Modified function outputs
  n=0;
  vector<MX> mfcn_out;  
  mfcn_out.push_back(g_z);                             mod_g_ = n++;
  
  // Constraint function
  MXFunction gfcn(mfcn_in,mfcn_out);
  gfcn.init();

  // Jacobian of the constraints
  MXFunction jac_fcn(mfcn_in,gfcn.jac(mod_x_,mod_g_));
  jac_fcn.init();
  log("Formed Jacobian of the constraints.");
  
  // Generate c code and load as DLL
  if(codegen_){
    dynamicCompilation(jac_fcn,jac_fcn_,"jac_fcn","Jacobian function");
  } else {
    jac_fcn_ = jac_fcn;
  }

  // Add multipliers to function inputs
  if(!gauss_newton_){
    n = mfcn_in.size();
    mfcn_in.push_back(g_lam);                          mod_g_lam_ = n++;
    mfcn_in.push_back(x_lam);                          mod_x_lam_ = n++;
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      mfcn_in.push_back(it->d_lam);                    it->mod_lam = n++;
    }
  }

  // Add gradient of the Lagrangian
  n = mfcn_out.size();
  mfcn_out.push_back(obj_z);                           mod_obj_ = n++;
  mfcn_out.push_back(gL_z);                            mod_gl_ = n++;  
  
  // Lagrangian gradient function
  MXFunction lgrad(mfcn_in,mfcn_out);
  lgrad.init();
  
  // Hessian of the Lagrangian
  MXFunction hes_fcn(mfcn_in,lgrad.jac(mod_x_,mod_gl_,false,!gauss_newton_));
  hes_fcn.init();
  log("Formed Hessian of the Lagrangian.");

  if(codegen_){
    dynamicCompilation(hes_fcn,hes_fcn_,"hes_fcn","Hessian function");
  } else {
    hes_fcn_ = hes_fcn;
  }

  // Definition of intermediate variables
  n = mfcn_out.size();
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    mfcn_out.push_back(it->d_def);         it->mod_def = n++;
    if(!gauss_newton_){
      mfcn_out.push_back(it->d_defL);      it->mod_defL = n++;
    }
  }

  // Modifier function
  MXFunction mfcn(mfcn_in,mfcn_out);
  mfcn.init();

  // Directional derivative of Z
  vector<vector<MX> > mfcn_fwdSeed(1,mfcn_in), mfcn_fwdSens(1,mfcn_out);
  vector<vector<MX> > mfcn_adjSeed,            mfcn_adjSens;

  // Linearization in the d-direction (see Equation (2.12) in Alberspeyer2010)
  fill(mfcn_fwdSeed[0].begin(),mfcn_fwdSeed[0].end(),MX());
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    mfcn_fwdSeed[0][it->mod_var] = it->d;
    if(!gauss_newton_){
      mfcn_fwdSeed[0][it->mod_lam] = it->d_lam;
    }
  }
  mfcn.eval(mfcn_in,mfcn_out,mfcn_fwdSeed,mfcn_fwdSens,mfcn_adjSeed,mfcn_adjSens,true);
  
  // Vector(s) b in Lifted Newton
  MX b_obj = mfcn_fwdSens[0][mod_gl_];
  MX b_g = mfcn_fwdSens[0][mod_g_];
  
  // Make sure that the vectors are dense
  makeDense(b_obj);
  makeDense(b_g);
  
  // Tangent function
  vector<MX> tan_fcn_out;
  n=0;
  tan_fcn_out.push_back(b_obj);                             tan_b_obj_ = n++;
  tan_fcn_out.push_back(b_g);                               tan_b_g_ = n++;  
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    tan_fcn_out.push_back(mfcn_fwdSens[0][it->mod_def]);    it->tan_lin = n++;
  }

  if(!gauss_newton_){
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      tan_fcn_out.push_back(mfcn_fwdSens[0][it->mod_defL]); it->tan_linL = n++;
    }
  }
  casadi_assert(n==tan_fcn_out.size());
  
  MXFunction tan_fcn(mfcn_in,tan_fcn_out);
  tan_fcn.setOption("name","tan_fcn");
  tan_fcn.setOption("number_of_fwd_dir",0);
  tan_fcn.setOption("number_of_adj_dir",0);
  tan_fcn.init();
  if(verbose_){
    cout << "Generated linearization function ( " << tan_fcn.getAlgorithmSize() << " nodes)." << endl;
  }
  
  // Generate c code and load as DLL
  if(codegen_){
    dynamicCompilation(tan_fcn,tan_fcn_,"tan_fcn","linearization function");
  } else {
    tan_fcn_ = tan_fcn;
  }

  // Expression a + A*du in Lifted Newton (Section 2.1 in Alberspeyer2010)
  MX du = msym("du",nx_);   // Step in u
  MX g_dlam;               // Step lambda_g
  if(!gauss_newton_){
    g_dlam = msym("g_dlam",g_lam.sparsity());
  }
  
  // Interpret the Jacobian-vector multiplication as a forward directional derivative
  fill(mfcn_fwdSeed[0].begin(),mfcn_fwdSeed[0].end(),MX());
  mfcn_fwdSeed[0][mod_x_] = du;
  if(!gauss_newton_){
    mfcn_fwdSeed[0][mod_g_lam_] = g_dlam;
  }
  mfcn.eval(mfcn_in,mfcn_out,mfcn_fwdSeed,mfcn_fwdSens,mfcn_adjSeed,mfcn_adjSens,true);    
  
  // Step expansion function inputs
  n = mfcn_in.size();
  mfcn_in.push_back(du);                                 mod_du_ = n++;
  if(!gauss_newton_){
    mfcn_in.push_back(g_dlam);                           mod_dlam_g_ = n++;
  }
    
  // Step expansion function outputs
  vector<MX> exp_fcn_out;
  n=0;
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    exp_fcn_out.push_back(mfcn_fwdSens[0][it->mod_def]); it->exp_def = n++;
  }

  if(!gauss_newton_){
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      exp_fcn_out.push_back(mfcn_fwdSens[0][it->mod_defL]); it->exp_defL = n++;
    }
  }
  
  // Step expansion function
  MXFunction exp_fcn(mfcn_in,exp_fcn_out);
  exp_fcn.setOption("number_of_fwd_dir",0);
  exp_fcn.setOption("number_of_adj_dir",0);
  exp_fcn.setOption("name","exp_fcn");
  exp_fcn.init();
  if(verbose_){
    cout << "Generated step expansion function ( " << exp_fcn.getAlgorithmSize() << " nodes)." << endl;
  }
  
  // Generate c code and load as DLL
  if(codegen_){
    dynamicCompilation(exp_fcn,exp_fcn_,"exp_fcn","step expansion function");
  } else {
    exp_fcn_ = exp_fcn;
  }  
  
  // Allocate QP data
  CRSSparsity sp_tr_B_obj = hes_fcn_.output().sparsity().transpose();
  qpH_ = DMatrix(sp_tr_B_obj.patternProduct(sp_tr_B_obj));
  qpA_ = jac_fcn_.output();
  qpG_.resize(nx_);
  qpB_.resize(ng_);

  // Allocate a QP solver
  QPSolverCreator qp_solver_creator = getOption("qp_solver");
  qp_solver_ = qp_solver_creator(qpH_.sparsity(),qpA_.sparsity());
  
  // Set options if provided
  if(hasSetOption("qp_solver_options")){
    Dictionary qp_solver_options = getOption("qp_solver_options");
    qp_solver_.setOption(qp_solver_options);
  }
  
  // Initialize the QP solver
  qp_solver_.init();
  if(verbose_){
    cout << "Allocated QP solver." << endl;
  }
  
  // Residual
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    it->res.resize(it->d.size(),0);
  }

  if(!gauss_newton_){
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      it->resL.resize(it->d_lam.size(),0);
    }
  }
  
  if(verbose_){
    cout << "NLP preparation completed" << endl;
  }
  
  // Header
  if(bool(getOption("print_header"))){
    cout << "-------------------------------------------" << endl;
    cout << "This is CasADi::SCPgen." << endl;
    if(gauss_newton_) {
      cout << "Using Gauss-Newton Hessian" << endl;
    } else {
      cout << "Using exact Hessian" << endl;
    }

    // Count the total number of variables
    int n_lifted = 0;
    for(vector<Var>::const_iterator i=v_.begin(); i!=v_.end(); ++i){
      n_lifted += i->n;
    }

    cout << endl;
    cout << "Number of reduced variables:               " << setw(9) << nx_ << endl;
    cout << "Number of reduced constraints:             " << setw(9) << ng_ << endl;
    cout << "Number of lifted variables/constraints:    " << setw(9) << n_lifted << endl;
    cout << "Number of parameters:                      " << setw(9) << np_ << endl;
    cout << "Total number of variables:                 " << setw(9) << (nx_+n_lifted) << endl;
    cout << "Total number of constraints:               " << setw(9) << (ng_+n_lifted) << endl;
    cout << endl;
    cout << "Iteration options:" << endl;

    cout << "{ \"maxiter\":" << maxiter_ << ", ";
    cout << "\"maxiter_ls\":" << maxiter_ls_ << ", ";
    cout << "\"c1\":" << c1_ << ", ";
    cout << "\"beta\":" << beta_ << ", ";
    cout << "\"merit_memsize\":" << merit_memsize_ << ", ";
    cout << "\"regularize\":" << regularize_ << ", ";
    cout << endl << "  ";
    cout << "\"tol_pr\":" << tol_pr_ << ", ";
    cout << "\"tol_du\":" << tol_du_ << ", ";
    cout << "\"tol_reg\":" << tol_reg_ << ", ";
    cout << "\"reg_threshold\":" << reg_threshold_ << "}" << endl;
    cout << endl;
  }
}

void SCPgenInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);

  checkInitialBounds();
  
  // Get problem data
  const vector<double>& x_init = input(NLP_X_INIT).data();
  const vector<double>& lbx = input(NLP_LBX).data();
  const vector<double>& ubx = input(NLP_UBX).data();
  const vector<double>& lbg = input(NLP_LBG).data();
  const vector<double>& ubg = input(NLP_UBG).data();  
  
  copy(x_init.begin(),x_init.end(),x_init_.begin());
  copy(lbx.begin(),lbx.end(),x_lb_.begin());
  copy(ubx.begin(),ubx.end(),x_ub_.begin());
  copy(lbg.begin(),lbg.end(),g_lb_.begin());
  copy(ubg.begin(),ubg.end(),g_ub_.begin());
  
  if(v_.size()>0){
    // Initialize lifted variables using the generated function
    vinit_fcn_.setInput(x_init_,0);
    vinit_fcn_.setInput(input(NLP_P),1);
    vinit_fcn_.evaluate();    
    for(int i=0; i<v_.size(); ++i){
      vinit_fcn_.getOutput(v_[i].init,i);
    }
  }
  if(verbose_){
    cout << "Passed initial guess" << endl;
  }

  // Reset dual guess
  if(!gauss_newton_){
    fill(g_lam_.begin(),g_lam_.end(),0);
    fill(g_dlam_.begin(),g_dlam_.end(),0);    
    fill(x_lam_.begin(),x_lam_.end(),0);
    fill(x_dlam_.begin(),x_dlam_.end(),0);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      fill(it->lam.begin(),it->lam.end(),0);
      fill(it->dlam.begin(),it->dlam.end(),0);
    }
  }
  
  // Objective value
  obj_k_ = numeric_limits<double>::quiet_NaN();

  // Reset line-search
  fill(merit_mem_.begin(),merit_mem_.end(),0.0);
  merit_ind_ = 0;

  // Current guess for the primal solution
  copy(x_init_.begin(),x_init_.end(),x_opt_.begin());
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    copy(it->init.begin(),it->init.end(),it->opt.begin());
  }

  // Get current time and reset timers
  double time1 = clock();
  t_eval_jac_ = t_eval_hes_ = t_eval_res_ = t_eval_tan_ = t_eval_exp_ = t_solve_qp_ = 0;

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
  
  // Reset iteration message
  iteration_note_ = string();

  // MAIN OPTIMZATION LOOP
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
    bool converged = pr_inf <= tol_pr_ && pr_step_ <= tol_pr_step_ && reg_ <= tol_reg_;
    if(gauss_newton_){
      converged = converged && iter!=0;
    } else {
      converged = converged && du_inf <= tol_du_;
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
    eval_tan();

    // Evaluate the constraint Jacobian
    eval_jac();
    
    // Evaluate the condensed Hessian
    eval_hess();
    
    // Regularize the QP
    if(regularize_){
      regularize();
    }

    // Solve the condensed QP
    solve_qp();

    // Expand the step
    eval_exp();
  
    // Line-search to take the step
    line_search(ls_iter, ls_success);
  }

  double time2 = clock();
  t_mainloop_ = double(time2-time1)/CLOCKS_PER_SEC;
  
  // Store optimal value
  cout << "optimal cost = " << obj_k_ << endl;

  // Save results to outputs
  output(NLP_COST).set(obj_k_);
  output(NLP_X_OPT).set(x_opt_);
  if(!gauss_newton_){
    output(NLP_LAMBDA_G).set(g_lam_);
    output(NLP_LAMBDA_X).set(x_lam_);
  }
  output(NLP_G).set(g_);
  
  // Write timers
  if(print_time_){
    cout << endl;
    cout << "time spent in eval_hes:    " << setw(9) << t_eval_hes_ << " s." << endl;
    cout << "time spent in eval_jac:    " << setw(9) << t_eval_jac_ << " s." << endl;
    cout << "time spent in eval_res:    " << setw(9) << t_eval_res_ << " s." << endl;
    cout << "time spent in eval_tan:    " << setw(9) << t_eval_tan_ << " s." << endl;
    cout << "time spent in eval_exp:    " << setw(9) << t_eval_exp_ << " s." << endl;
    cout << "time spent in solve_qp:    " << setw(9) << t_solve_qp_ << " s." << endl;
    cout << "time spent in main loop:   " << setw(9) << t_mainloop_ << " s." << endl;
  }

  // Save statistics
  stats_["iter_count"] = iter;

  cout << endl;
}  

void SCPgenInternal::dynamicCompilation(FX& f, FX& f_gen, std::string fname, std::string fdescr){
#ifdef WITH_DL 

  // C compiler command
  string compiler = getOption("compiler");

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
  string compile_command = compiler + " " + dlflag + " " + cname + " -o " + dlname;
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

#else // WITH_DL 
    casadi_error("Codegen in SCPgen requires CasADi to be compiled with option \"WITH_DL\" enabled");
#endif // WITH_DL 
}

double SCPgenInternal::primalInfeasibility(){
  // L1-norm of the primal infeasibility
  double pr_inf = 0;
  
  // Simple bounds
  for(int i=0; i<nx_; ++i) pr_inf +=  ::fmax(x_opt_[i]-x_ub_[i],0.);
  for(int i=0; i<nx_; ++i) pr_inf +=  ::fmax(x_lb_[i]-x_opt_[i],0.);

  // Lifted variables
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    for(int i=0; i<it->n; ++i) pr_inf += ::fabs(it->res[i]);
  }
  
  // Nonlinear bounds
  for(int i=0; i<ng_; ++i) pr_inf += ::fmax(g_[i]-g_ub_[i],0.);
  for(int i=0; i<ng_; ++i) pr_inf += ::fmax(g_lb_[i]-g_[i],0.);
  
  return pr_inf;
}  

double SCPgenInternal::dualInfeasibility(){
  // Not implemented for Gauss-Newton
  if(gauss_newton_) return 0;

  // L1-norm of the dual infeasibility
  double du_inf = 0;
  
  // Lifted variables
  for(int i=0; i<ngL_; ++i) du_inf += ::fabs(gL_[i]);

  return du_inf;
}

void SCPgenInternal::printIteration(std::ostream &stream){
  stream << setw(4)  << "iter";
  stream << setw(14) << "objective";
  stream << setw(11) << "inf_pr";
  stream << setw(11) << "inf_du";
  stream << setw(11) << "pr_step";
  stream << setw(11) << "du_step";
  stream << setw(8) << "lg(rg)";
  stream << setw(3) << "ls";
  stream << ' ';

  // Print variables
  for(vector<int>::const_iterator i=print_x_.begin(); i!=print_x_.end(); ++i){
    stream << setw(9) << name_x_.at(*i);
  }
  
  stream << endl;
  stream.unsetf( std::ios::floatfield);
}
  
void SCPgenInternal::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double rg, int ls_trials, bool ls_success){
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

  // Print variables
  for(vector<int>::const_iterator i=print_x_.begin(); i!=print_x_.end(); ++i){
    stream << setw(9) << setprecision(4) << x_opt_.at(*i);
  }

  // Print note
  if(!iteration_note_.empty()){
    stream << "   " << iteration_note_;
    iteration_note_ = string();
  }

  stream.unsetf( std::ios::floatfield);
  stream << endl;
}

void SCPgenInternal::eval_jac(){
  // Get current time
  double time1 = clock();

  // Pass current parameter guess
  jac_fcn_.setInput(input(NLP_P),mod_p_);

  // Pass primal step/variables
  jac_fcn_.setInput(x_opt_, mod_x_);
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    jac_fcn_.setInput(it->res, it->mod_var);
  }

  // Evaluate condensed Jacobian
  jac_fcn_.evaluate();

  // Get the Jacobian
  jac_fcn_.getOutput(qpA_);

  double time2 = clock();
  t_eval_jac_ += double(time2-time1)/CLOCKS_PER_SEC;
}

void SCPgenInternal::eval_hess(){
  // Get current time
  double time1 = clock();

  // Pass parameters
  hes_fcn_.setInput(input(NLP_P),mod_p_);

  // Pass primal step/variables  
  hes_fcn_.setInput(x_opt_, mod_x_);
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    hes_fcn_.setInput(it->res, it->mod_var);
  }

  // Pass dual steps/variables
  if(!gauss_newton_){
    hes_fcn_.setInput(g_lam_,mod_g_lam_);
    hes_fcn_.setInput(x_lam_, mod_x_lam_);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      hes_fcn_.setInput(it->resL, it->mod_lam);
    }
  }

  // Evaluate condensed Hessian
  hes_fcn_.evaluate();

  // Get the objective function terms in the QP
  if(gauss_newton_){
    // Hessian of the lagrangian
    const DMatrix& B_obj =  hes_fcn_.output();
    fill(qpH_.begin(),qpH_.end(),0);
    DMatrix::mul_no_alloc_tn(B_obj,B_obj,qpH_);

    // Gradient of the objective
    fill(qpG_.begin(),qpG_.end(),0);
    DMatrix::mul_no_alloc_tn(B_obj,gL_,qpG_);
  } else {
    // Hessian of the lagrangian
    hes_fcn_.getOutput(qpH_);

    // Gradient of the lagrangian
    copy(gL_.begin(),gL_.end(),qpG_.begin());

    // Remove the contribution from the simple bounds multipliers
    for(int i=0; i<nx_; ++i){
      qpG_[i] -= x_lam_[i];
    }

    // Remove the contribution from the nonlinear multipliers to get the gradient of the objective
    const vector<double> &qpA_data = qpA_.data();
    const vector<int> &qpA_rowind = qpA_.rowind();
    const vector<int> &qpA_col = qpA_.col();
    for(int i=0; i<ng_; ++i){
      for(int el=qpA_rowind[i]; el<qpA_rowind[i+1]; ++el){
	int j=qpA_col[el];
	qpG_[j] -= qpA_data[el]*g_lam_[i];
      }
    }
  }

  double time2 = clock();
  t_eval_hes_ += double(time2-time1)/CLOCKS_PER_SEC;
}

void SCPgenInternal::eval_res(){
  // Get current time
  double time1 = clock();

  // Pass parameters
  res_fcn_.setInput(input(NLP_P),res_p_);

  // Pass primal variables to the residual function for initial evaluation
  res_fcn_.setInput(x_opt_,res_x_);
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    res_fcn_.setInput(it->opt,it->res_var);
  }

  // Pass dual variables to the residual function for initial evaluation
  if(!gauss_newton_){
    res_fcn_.setInput(g_lam_,res_g_lam_);
    res_fcn_.setInput(x_lam_,res_x_lam_);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      res_fcn_.setInput(it->lam,it->res_lam);
    }
  }
  
  // Evaluate residual function
  res_fcn_.evaluate();

  // Get objective
  obj_k_ = res_fcn_.output(res_obj_).toScalar();

  // Get objective gradient
  res_fcn_.getOutput(gL_,res_gl_);

  // Get constraints
  res_fcn_.getOutput(g_,res_g_);

  // Get residuals
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    res_fcn_.getOutput(it->res,  it->res_d);
    if(!gauss_newton_){
      res_fcn_.getOutput(it->resL, it->res_lam_d);
    }
  }
  
  // Parameter sensitivities
  res_fcn_.getOutput(output(NLP_LAMBDA_P),res_p_d_);

  double time2 = clock();
  t_eval_res_ += double(time2-time1)/CLOCKS_PER_SEC;
}

void SCPgenInternal::eval_tan(){
  // Get current time
  double time1 = clock();

  // Pass current parameter guess
  tan_fcn_.setInput(input(NLP_P),mod_p_);

  // Pass primal step/variables
  tan_fcn_.setInput(x_opt_, mod_x_);
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    tan_fcn_.setInput(it->res, it->mod_var);
  }

  // Pass dual steps/variables
  if(!gauss_newton_){
    tan_fcn_.setInput(g_lam_,mod_g_lam_);
    tan_fcn_.setInput(x_lam_, mod_x_lam_);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      tan_fcn_.setInput(it->resL, it->mod_lam);
    }
  }

  // Evaluate to get QP
  tan_fcn_.evaluate();

  // Get condensed vectors
  transform(g_.begin(),g_.end(),tan_fcn_.output(tan_b_g_).begin(),qpB_.begin(),std::minus<double>());
  transform(gL_.begin(),gL_.end(),tan_fcn_.output(tan_b_obj_).begin(),gL_.begin(),std::minus<double>());
  
  // Expanded primal step (first part)
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    const DMatrix& dv = tan_fcn_.output(it->tan_lin);
    copy(dv.begin(),dv.end(),it->step.begin());
  }
  
  // Expanded dual step (first part)
  if(!gauss_newton_){
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      const DMatrix& dlam_v = tan_fcn_.output(it->tan_linL);
      copy(dlam_v.begin(),dlam_v.end(),it->dlam.begin());
    }
  }

  double time2 = clock();
  t_eval_tan_ += double(time2-time1)/CLOCKS_PER_SEC;
}

void SCPgenInternal::regularize(){
  casadi_assert(nx_==2);
  
  // Regularization
  reg_ = 0;
  
  // Check the smallest eigenvalue of the Hessian
  double a = qpH_.elem(0,0);
  double b = qpH_.elem(0,1);
  double c = qpH_.elem(1,0);
  double d = qpH_.elem(1,1);
  
  // Make sure no not a numbers
  casadi_assert(a==a && b==b && c==c &&  d==d);
  
  // Make sure symmetric
  if(b!=c){
    casadi_assert_warning(fabs(b-c)<1e-10,"Hessian is not symmetric: " << b << " != " << c);
    qpH_.elem(1,0) = c = b;
  }
  
  double eig_smallest = (a+d)/2 - std::sqrt(4*b*c + (a-d)*(a-d))/2;
  if(eig_smallest<reg_threshold_){
    // Regularization
    reg_ = reg_threshold_-eig_smallest;
    qpH_(0,0) += reg_;
    qpH_(1,1) += reg_;
  }
}

void SCPgenInternal::solve_qp(){  
  // Get current time
  double time1 = clock();

  // Solve the QP
  qp_solver_.setInput(qpH_,QP_H);
  qp_solver_.setInput(qpG_,QP_G);
  qp_solver_.setInput(qpA_,QP_A);
  std::transform(x_lb_.begin(),x_lb_.end(), x_opt_.begin(),qp_solver_.input(QP_LBX).begin(),std::minus<double>());
  std::transform(x_ub_.begin(),x_ub_.end(), x_opt_.begin(),qp_solver_.input(QP_UBX).begin(),std::minus<double>()); 
  std::transform(g_lb_.begin(),g_lb_.end(), qpB_.begin(),qp_solver_.input(QP_LBA).begin(),std::minus<double>());
  std::transform(g_ub_.begin(),g_ub_.end(), qpB_.begin(),qp_solver_.input(QP_UBA).begin(),std::minus<double>());
  
  qp_solver_.evaluate();
  
  // Condensed primal step
  const DMatrix& du = qp_solver_.output(QP_PRIMAL);
  copy(du.begin(),du.end(),x_step_.begin());
  
  if(!gauss_newton_){

    // Condensed dual step (simple bounds)
    const DMatrix& lam_x_new = qp_solver_.output(QP_LAMBDA_X);
    copy(lam_x_new.begin(),lam_x_new.end(),x_dlam_.begin());
    std::transform(x_dlam_.begin(),x_dlam_.end(),x_lam_.begin(),x_dlam_.begin(),std::minus<double>());

    // Condensed dual step (nonlinear bounds)
    const DMatrix& lam_g_new = qp_solver_.output(QP_LAMBDA_A);
    copy(lam_g_new.begin(),lam_g_new.end(),g_dlam_.begin());
    std::transform(g_dlam_.begin(),g_dlam_.end(),g_lam_.begin(),g_dlam_.begin(),std::minus<double>());
  }

  double time2 = clock();
  t_solve_qp_ += double(time2-time1)/CLOCKS_PER_SEC;
}

void SCPgenInternal::line_search(int& ls_iter, bool& ls_success){
  // Make sure that we have a decent direction 
  if(!gauss_newton_){
    // Get the curvature in the step direction
    double gain = DMatrix::quad_form(qpH_,x_step_);
    if (gain < 0){
      iteration_note_ = "Hessian indefinite in the search direction";
    }
  }

  // Calculate penalty parameter of merit function
  sigma_ = 0;
  sigma_ = std::max(sigma_,1.01*norm_inf(qp_solver_.output(QP_LAMBDA_X).data()));
  sigma_ = std::max(sigma_,1.01*norm_inf(qp_solver_.output(QP_LAMBDA_A).data()));
  
  // Calculate L1-merit function in the actual iterate
  double l1_infeas = primalInfeasibility();

  // Right-hand side of Armijo condition
  double F_sens = 0;
  for(int i=0; i<nx_; ++i) F_sens += x_step_[i] * qpG_[i];
  double L1dir = F_sens - sigma_ * l1_infeas;
  double L1merit = obj_k_ + sigma_ * l1_infeas;
  
  // Storing the actual merit function value in a list
  merit_mem_[merit_ind_] = L1merit;
  ++merit_ind_ %= merit_memsize_;

  // Stepsize
  double t = 1.0, t_prev = 0.0;
  double fk_cand;
  
  // Merit function value in candidate
  double L1merit_cand = 0;
  
  // Reset line-search counter, success marker
  ls_iter = 0;
  ls_success = false;
  
  // Line-search
  //log("Starting line-search");
  
  // Line-search loop
  while (true){
    
    // Take the primal step
    for(int i=0; i<nx_; ++i) x_opt_[i] += (t-t_prev) * x_step_[i];
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      for(int i=0; i<it->n; ++i) it->opt[i] += (t-t_prev) * it->step[i];
    }
    
    // Take the dual step
    if(!gauss_newton_){
      for(int i=0; i<ng_; ++i) 	g_lam_[i] += (t-t_prev) * g_dlam_[i];
      for(int i=0; i<nx_; ++i)  x_lam_[i] += (t-t_prev) * x_dlam_[i];
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
	for(int i=0; i<it->n; ++i) it->lam[i] += (t-t_prev) * it->dlam[i];
      }
    }
    
    // Evaluate residual function to get objective and constraints (and residuals for the next iteration)
    eval_res();
    ls_iter++;      
    
    // Calculating merit-function in candidate
    l1_infeas = primalInfeasibility();
    L1merit_cand = obj_k_ + sigma_ * l1_infeas;
    
    // Calculating maximal merit function value so far
    double meritmax = *max_element(merit_mem_.begin(), merit_mem_.end());
    if (L1merit_cand <= meritmax + t * c1_ * L1dir){
      
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
  for(vector<double>::const_iterator i=x_step_.begin(); i!=x_step_.end(); ++i) pr_step_ += fabs(*i);
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    for(vector<double>::const_iterator i=it->step.begin(); i!=it->step.end(); ++i) pr_step_ += fabs(*i);
  }
  pr_step_ *= t;

  // Calculate the dual step-size
  if(!gauss_newton_){
    du_step_ = 0;
    for(vector<double>::const_iterator i=g_dlam_.begin(); i!=g_dlam_.end(); ++i){
      du_step_ += fabs(*i);
    }

    for(vector<double>::const_iterator i=x_dlam_.begin(); i!=x_dlam_.end(); ++i) du_step_ += fabs(*i);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      for(vector<double>::const_iterator i=it->dlam.begin(); i!=it->dlam.end(); ++i) du_step_ += fabs(*i);
    }
    du_step_ *= t;
  }
}

void SCPgenInternal::eval_exp(){
  // Get current time
  double time1 = clock();

  // Pass current parameter guess
  exp_fcn_.setInput(input(NLP_P), mod_p_);

  // Pass primal step/variables
  exp_fcn_.setInput(x_step_, mod_du_);
  exp_fcn_.setInput(x_opt_,mod_x_);
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    exp_fcn_.setInput(it->res,it->mod_var);
  }
  
  // Pass dual step/variables
  if(!gauss_newton_){
    exp_fcn_.setInput(g_dlam_,mod_dlam_g_);
    exp_fcn_.setInput(g_lam_,mod_g_lam_);
    exp_fcn_.setInput(x_lam_,mod_x_lam_);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      exp_fcn_.setInput(it->resL,it->mod_lam);
    }
  }

  // Perform the step expansion
  exp_fcn_.evaluate();

  // Expanded primal step (second part)
  for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
    const DMatrix& dv = exp_fcn_.output(it->exp_def);
    transform(dv.begin(),dv.end(),it->step.begin(),it->step.begin(),std::minus<double>());
  }
  
  // Expanded dual step (second part)
  if(!gauss_newton_){
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      const DMatrix& dlam_v = exp_fcn_.output(it->exp_defL);
      transform(dlam_v.begin(),dlam_v.end(),it->dlam.begin(),it->dlam.begin(),std::minus<double>());
    }
  }

  double time2 = clock();
  t_eval_exp_ += double(time2-time1)/CLOCKS_PER_SEC;
}
  

} // namespace CasADi
