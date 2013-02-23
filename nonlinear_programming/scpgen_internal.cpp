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
  addOption("c1",                OT_REAL,       1E-4,             "Armijo condition, coefficient of decrease in merit");
  addOption("beta",              OT_REAL,       0.8,              "Line-search parameter, restoration factor of stepsize");
  addOption("merit_memory",      OT_INTEGER,      4,              "Size of memory to store history of merit function values");
  addOption("lbfgs_memory",      OT_INTEGER,     10,              "Size of L-BFGS memory.");
  addOption("regularize",        OT_BOOLEAN,  false,              "Automatic regularization of Lagrange Hessian.");
  addOption("print_header",      OT_BOOLEAN,   true,              "Print the header with problem statistics");
  addOption("codegen",           OT_BOOLEAN,  false,              "C-code generation");
  addOption("reg_threshold",     OT_REAL,      1e-8,              "Threshold for the regularization.");
  addOption("name_x",      OT_STRINGVECTOR,  GenericType(),       "Names of the variables.");
  addOption("print_x",           OT_INTEGERVECTOR,  GenericType(), "Which variables to print.");
  
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
  merit_memsize_ = getOption("merit_memory");
  lbfgs_memory_ = getOption("lbfgs_memory");
  tol_pr_ = getOption("tol_pr");
  tol_du_ = getOption("tol_du");
  regularize_ = getOption("regularize");
  codegen_ = getOption("codegen");
  reg_threshold_ = getOption("reg_threshold");

  // Name the components
  if(hasSetOption("name_x")){
    name_x_ = getOption("name_x");
    casadi_assert(name_x_.size()==n_);
  } else {
    stringstream ss;
    name_x_.resize(n_);
    for(int i=0; i<n_; ++i){
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
  vector<MX> x = vdef_fcn.inputExpr();
  vector<MX> f = vdef_fcn.outputExpr();

  // Get the dimensions
  x_.resize(x.size());
  for(int i=0; i<x_.size(); ++i){
    x_[i].v = x[i];
    x_[i].n = x[i].size();
  }

  // Allocate memory
  lbu_.resize(n_,-numeric_limits<double>::infinity());
  ubu_.resize(n_, numeric_limits<double>::infinity());
  lbg_.resize(m_,-numeric_limits<double>::infinity());
  ubg_.resize(m_, numeric_limits<double>::infinity());
  g_.resize(m_, numeric_limits<double>::quiet_NaN());
  if(!gauss_newton_){
    lambda_g_.resize(m_,0);
    dlambda_g_.resize(m_,0);
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
    obj = inner_prod(f[0],f[0])/2;
    lam_con[0] = f[0];
    ngL_ = lam_con[0].size();

  } else {
    // Scalar objective function
    obj = f[0];
    
    // Lagrange multipliers for the simple bounds on u and Lagrange multipliers corresponding to the definition of the dependent variables
    stringstream ss;
    for(int i=0; i<x_.size(); ++i){
      ss.str(string());
      ss << "lam_x" << i;
      lam[i] = msym(ss.str(),x[i].sparsity());
    }

    // Lagrange multipliers for the nonlinear constraints
    lam_g = msym("lam_g",m_);

    if(verbose_){
      cout << "Allocated intermediate variables." << endl;
    }
   
    // Adjoint sweep to get the definitions of the lifted dual variables (Equation 3.8 in Albersmeyer2010)
    vector<vector<MX> > fseed,fsens,aseed(1),asens(1);
    aseed[0].push_back(1.0);
    aseed[0].push_back(lam_g);
    aseed[0].insert(aseed[0].end(),lam.begin()+2,lam.end());
    vdef_fcn.eval(x,f,fseed,fsens,aseed,asens,true);

    for(int i=0; i<x_.size(); ++i){
      lam_con[i] = asens[0].at(i);
      if(lam_con[i].isNull()){
	lam_con[i] = MX(x[i].sparsity());
      }
    }
    ngL_ = n_;

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
    res_fcn_in.push_back(x[i]);        x_[i].res_var = n++;
    if(!gauss_newton_){
      res_fcn_in.push_back(lam[i]);    x_[i].res_lam = n++;
    }
  }

  // Outputs
  vector<MX> res_fcn_out;
  n=0;
  res_fcn_out.push_back(obj);                         res_obj_ = n++;
  res_fcn_out.push_back(lam_con[0]);                  res_gl_ = n++;
  res_fcn_out.push_back(f[1]);                        res_g_ = n++;
  for(int i=1; i<x_.size(); ++i){
    res_fcn_out.push_back(i==1 ? -x[i] : f[i]-x[i]);  x_[i].res_d = n++;

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
  if(codegen_){
    dynamicCompilation(res_fcn,res_fcn_,"res_fcn","residual function");
  } else {
    res_fcn_ = res_fcn;
  }

  // Declare difference vector d and substitute out p and v
  vector<MX> d;
  vector<MX> d_def;
  stringstream ss;
  for(int i=1; i<x_.size(); ++i){
    ss.str(string());
    ss << "d" << i;
    MX d_i = msym(ss.str(),x[i].sparsity());
    d.push_back(d_i);    
    d_def.push_back(i==1 ? x[i]-d_i : f[i]-d_i);
  }

  // Declare difference vector lam_d and substitute out lam
  vector<MX> lam_d;
  vector<MX> lam_d_def;
  if(!gauss_newton_){
    for(int i=1; i<x_.size(); ++i){
      ss.str(string());
      ss << "lam_d" << i;
      MX lam_d_i = msym(ss.str(),x[i].sparsity());
      lam_d.push_back(lam_d_i);
      lam_d_def.push_back(lam_con[i]-lam_d_i);
    }
  }

  // Variables to be substituted and their definitions
  vector<MX> svar, sdef;
  svar.insert(svar.end(),x.begin()+1,x.end());
  sdef.insert(sdef.end(),d_def.begin(),d_def.end());
  if(!gauss_newton_){
    svar.insert(svar.end(),lam.rbegin(),lam.rend()-1);
    sdef.insert(sdef.end(),lam_d_def.rbegin(),lam_d_def.rend());    
  }

  vector<MX> ex(3);
  ex[0] = obj;
  ex[1] = f[1];
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
  vector<MX> zfcn_in;
  n=0;
  zfcn_in.push_back(x[1]);                               z_p_ = n++;
  if(!gauss_newton_){
    zfcn_in.push_back(lam_g);                         z_lam_g_ = n++;
  }
  for(int i=0; i<x_.size(); ++i){
    zfcn_in.push_back(i==0 ? x[0] :     d[i-1]);    x_[i].z_var = n++;
    if(!gauss_newton_){
      zfcn_in.push_back(i==0 ? lam[0] : lam_d[i-1]);  x_[i].z_lam = n++;
    }
  }

  // Outputs
  n=0;
  vector<MX> zfcn_out;
  zfcn_out.push_back(obj_z);                   z_obj_ = n++;
  zfcn_out.push_back(gL_z);                    z_gl_ = n++;
  zfcn_out.push_back(g_z);                     z_g_ = n++;
  for(int i=1; i<x_.size(); ++i){
    zfcn_out.push_back(d_def[i-1]);            x_[i].z_def = n++;
    if(!gauss_newton_){
      zfcn_out.push_back(lam_d_def[i-1]);      x_[i].z_defL = n++;
    }
  }

  MXFunction zfcn(zfcn_in,zfcn_out);
  zfcn.setOption("name","zfcn");
  zfcn.init();
  if(verbose_){
    cout << "Generated reconstruction function ( " << zfcn.getAlgorithmSize() << " nodes)." << endl;
  }

  // Directional derivative of Z
  vector<vector<MX> > Z_fwdSeed(1,zfcn_in);
  vector<vector<MX> > Z_fwdSens(1,zfcn_out);
  vector<vector<MX> > Z_adjSeed;
  vector<vector<MX> > Z_adjSens;

  // Expression a + A*du in Lifted Newton (Section 2.1 in Alberspeyer2010)
  MX du = msym("du",n_);   // Step in u
  MX dlam_g;                // Step lambda_g
  if(!gauss_newton_){
    dlam_g = msym("dlam_g",lam_g.sparsity());
  }
  
  // Interpret the Jacobian-vector multiplication as a forward directional derivative
  fill(Z_fwdSeed[0].begin(),Z_fwdSeed[0].end(),MX());
  Z_fwdSeed[0][x_[0].z_var] = du;
  if(!gauss_newton_){
    Z_fwdSeed[0][z_lam_g_] = dlam_g;
  }
  zfcn.eval(zfcn_in,zfcn_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true);    
  
  // Step expansion function inputs
  vector<MX> exp_fcn_in;
  n=0;
  exp_fcn_in.push_back(x[1]);                            exp_p_ = n++;
  if(!gauss_newton_){
    exp_fcn_in.push_back(lam_g);                         exp_lam_g_ = n++;
  }
  exp_fcn_in.push_back(du);                              exp_du_ = n++;
  if(!gauss_newton_){
    exp_fcn_in.push_back(dlam_g);                        exp_dlam_g_ = n++;
  }
  
  for(int i=0; i<x_.size(); ++i){
    exp_fcn_in.push_back(   i==0 ? x[0] :     d[i-1]);   x_[i].exp_var = n++;
    if(!gauss_newton_){
      exp_fcn_in.push_back( i==0 ? lam[0] : lam_d[i-1]);   x_[i].exp_lam = n++;
    }
  }
  
  // Step expansion function outputs
  vector<MX> exp_fcn_out;
  n=0;
  exp_fcn_out.push_back(Z_fwdSens[0][z_obj_]);           exp_osens_ = n++;
  exp_fcn_out.push_back(Z_fwdSens[0][z_gl_]);            exp_curve_ = n++;
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
  if(codegen_){
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
  zfcn.eval(zfcn_in,zfcn_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true);   
  
  // Vector(s) b in Lifted Newton
  MX b_obj = Z_fwdSens[0][z_gl_];
  MX b_g = Z_fwdSens[0][z_g_];

  // Differentiate with respect to the step to get the matrix B in Lifted Newton
  MX B_obj = zfcn.jac(x_[0].z_var,z_gl_,false,!gauss_newton_); // Exploit Hessian symmetry
  if(verbose_){
    cout << "Formed B_obj (" << B_obj.size1() << "-by-" << B_obj.size2() << ", "<< B_obj.size() << " nnz)" << endl;
  }
  MX B_g = zfcn.jac(x_[0].z_var,z_g_);
  if(B_g.empty()){
    B_g = MX::sparse(0,n_);
  }
  if(verbose_){
    cout << "Formed B_g (" << B_g.size1() << "-by-" << B_g.size2() << ", "<< B_g.size() << " nnz)" << endl;
  }
  
  // Make sure that the vectors are dense
  makeDense(b_obj);
  makeDense(b_g);
  
  // Quadratic approximation
  vector<MX> qp_fcn_in;
  n=0;
  qp_fcn_in.push_back(x[1]);             qpf_p_ = n++;
  if(!gauss_newton_){
    qp_fcn_in.push_back(lam_g);          qpf_lam_g_ = n++;
  }
  for(int i=0; i<x_.size(); ++i){
    qp_fcn_in.push_back(i==0 ?     x[i] : d[i-1]);        x_[i].qpf_var = n++;
    if(!gauss_newton_){
      qp_fcn_in.push_back(i==0 ? lam[i] : lam_d[i-1]);    x_[i].qpf_lam = n++;
    }
  }
  casadi_assert(n==qp_fcn_in.size());  

  vector<MX> qp_fcn_out;
  n=0;
  qp_fcn_out.push_back(b_obj);       qpf_b_obj_ = n++;
  qp_fcn_out.push_back(B_obj);       qpf_B_obj_ = n++;
  qp_fcn_out.push_back(b_g);         qpf_b_g_ = n++;
  qp_fcn_out.push_back(B_g);         qpf_B_g_ = n++;
   for(int i=1; i<x_.size(); ++i){
    qp_fcn_out.push_back(Z_fwdSens[0][x_[i].z_def]);    x_[i].qpf_def = n++;
    if(!gauss_newton_){
      qp_fcn_out.push_back(Z_fwdSens[0][x_[i].z_defL]); x_[i].qpf_defL = n++;
    }
   }
  casadi_assert(n==qp_fcn_out.size());
  
  MXFunction qp_fcn(qp_fcn_in,qp_fcn_out);
  qp_fcn.setOption("name","qp_fcn");
  qp_fcn.setOption("number_of_fwd_dir",0);
  qp_fcn.setOption("number_of_adj_dir",0);
  qp_fcn.init();
  if(verbose_){
    cout << "Generated linearization function ( " << qp_fcn.getAlgorithmSize() << " nodes)." << endl;
  }
  
  // Generate c code and load as DLL
  if(codegen_){
    dynamicCompilation(qp_fcn,qp_fcn_,"qp_fcn","linearization function");
  } else {
    qp_fcn_ = qp_fcn;
  }
  
  // Allocate QP data
  CRSSparsity sp_tr_B_obj = B_obj.sparsity().transpose();
  qpH_ = DMatrix(sp_tr_B_obj.patternProduct(sp_tr_B_obj));
  qpA_ = DMatrix(B_g.sparsity());
  qpG_.resize(n_);
  qpB_.resize(m_);

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
  for(int i=1; i<x_.size(); ++i){
    x_[i].res.resize(d[i-1].size(),0);
    if(!gauss_newton_){
      x_[i].resL.resize(lam_d[i-1].size(),0);
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
    for(vector<Var>::const_iterator i=x_.begin()+2; i!=x_.end(); ++i){
      n_lifted += i->n;
    }

    cout << endl;
    cout << "Number of reduced variables:               " << setw(9) << n_ << endl;
    cout << "Number of reduced constraints:             " << setw(9) << m_ << endl;
    cout << "Number of lifted variables/constraints:    " << setw(9) << n_lifted << endl;
    cout << "Number of parameters:                      " << setw(9) << x_[1].n << endl;
    cout << "Total number of variables:                 " << setw(9) << (n_+n_lifted) << endl;
    cout << "Total number of constraints:               " << setw(9) << (m_+n_lifted) << endl;
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
  if(parametric_){
    const vector<double>& p = input(NLP_P).data();
    copy(p.begin(),p.end(),x_[1].init.begin());
  }

  copy(x_init.begin(),x_init.end(),x_[0].init.begin());
  copy(lbx.begin(),lbx.end(),lbu_.begin());
  copy(ubx.begin(),ubx.end(),ubu_.begin());
  copy(lbg.begin(),lbg.end(),lbg_.begin());
  copy(ubg.begin(),ubg.end(),ubg_.begin());
  
  if(x_.size()>2){
    // Initialize lifted variables using the generated function
    for(int i=0; i<2; ++i){
      vinit_fcn_.setInput(x_[i].init,i);
    }
    vinit_fcn_.evaluate();    
    for(int i=2; i<x_.size(); ++i){
      vinit_fcn_.getOutput(x_[i].init,i-2);
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
  
  double toldx_ = 1e-9;

  // Objective value
  obj_k_ = numeric_limits<double>::quiet_NaN();

  // Reset line-search
  merit_mem_.clear();

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
  
  // Store optimal value
  cout << "optimal cost = " << obj_k_ << endl;

  // Save results to outputs
  output(NLP_COST).set(obj_k_);
  output(NLP_X_OPT).set(x_[0].opt);
  if(!gauss_newton_){
    output(NLP_LAMBDA_G).set(lambda_g_);
    output(NLP_LAMBDA_X).set(x_[0].lam);
  }
  output(NLP_G).set(g_);
  
  // Save statistics
  stats_["iter_count"] = iter;
}  

void SCPgenInternal::dynamicCompilation(FX& f, FX& f_gen, std::string fname, std::string fdescr){
#ifdef WITH_DL 

  // C compiler
  string compiler = "clang";

  // Optimization flag
  //string oflag = "-O3";
  string oflag = "-O2";

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

#else // WITH_DL 
    casadi_error("Codegen in SCPgen requires CasADi to be compiled with option \"WITH_DL\" enabled");
#endif // WITH_DL 
}

double SCPgenInternal::primalInfeasibility(){
  // L1-norm of the primal infeasibility
  double pr_inf = 0;
  
  // Variable bounds
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    if(it==x_.begin()){
      
      // Simple bounds
      for(int i=0; i<it->n; ++i) pr_inf +=  ::fmax(it->opt[i]-ubu_[i],0.);
      for(int i=0; i<it->n; ++i) pr_inf +=  ::fmax(lbu_[i]-it->opt[i],0.);
    } else {
      
      // Parameters and other lifted variables
      for(int i=0; i<it->n; ++i) pr_inf += ::fabs(it->res[i]);
    }
  }
  
  // Nonlinear bounds
  for(int i=0; i<m_; ++i) pr_inf += ::fmax(g_[i]-ubg_[i],0.);
  for(int i=0; i<m_; ++i) pr_inf += ::fmax(lbg_[i]-g_[i],0.);
  
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
    stream << setw(9) << setprecision(4) << x_[0].opt.at(*i);
  }

  stream << scientific;
  stream << endl;
}

void SCPgenInternal::eval_res(){
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
  res_fcn_.getOutput(gL_,res_gl_);

  // Get constraints
  res_fcn_.getOutput(g_,res_g_);

  // Get residuals
  for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
    res_fcn_.getOutput(it->res,  it->res_d);
    if(!gauss_newton_){
      res_fcn_.getOutput(it->resL, it->res_lam_d);
    }
  }
  
  // Embedded parameters
  if(parametric_){
    const vector<double>& p = input(NLP_P).data();
    transform(x_[1].res.begin(),x_[1].res.end(),p.begin(),x_[1].res.begin(),std::plus<double>());
  }
}

void SCPgenInternal::eval_qpf(){
  // Pass current parameter guess
  qp_fcn_.setInput(x_[1].opt,qpf_p_);

  // Pass primal step/variables
  for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
    qp_fcn_.setInput(it==x_.begin() ? it->opt : it->res, it->qpf_var);
  }

  // Pass dual steps/variables
  if(!gauss_newton_){
    qp_fcn_.setInput(lambda_g_,qpf_lam_g_);
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      qp_fcn_.setInput(it==x_.begin() ? it->lam : it->resL, it->qpf_lam);
    }
  }

  // Evaluate to get QP
  qp_fcn_.evaluate();

  // Get condensed vectors
  transform(g_.begin(),g_.end(),qp_fcn_.output(qpf_b_g_).begin(),qpB_.begin(),std::minus<double>());
  transform(gL_.begin(),gL_.end(),qp_fcn_.output(qpf_b_obj_).begin(),gL_.begin(),std::minus<double>());
  
  // Get condensed linear constraint
  qp_fcn_.getOutput(qpA_,qpf_B_g_);
  
  // Get the reduced Hessian
  if(gauss_newton_){
    const DMatrix& B_obj =  qp_fcn_.output(qpf_B_obj_);
    fill(qpH_.begin(),qpH_.end(),0);
    DMatrix::mul_no_alloc_tn(B_obj,B_obj,qpH_);

    fill(qpG_.begin(),qpG_.end(),0);
    DMatrix::mul_no_alloc_tn(B_obj,gL_,qpG_);
  } else {
    qp_fcn_.getOutput(qpH_,qpf_B_obj_);
    copy(gL_.begin(),gL_.end(),qpG_.begin());
  }

  // Change the gradient of the Lagrangian to the gradient of the objective
  if(!gauss_newton_){
    // Subtract a multiple of the condensed Hessian to 
    for(int i=0; i<qpA_.size1(); ++i){
      for(int el=qpA_.rowind(i); el<qpA_.rowind(i+1); ++el){
	int j=qpA_.col(el);
	qpG_[j] -= qpA_.at(el)*lambda_g_[i];
      }
    }
  }

  // Expanded primal step (first part)
  for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
    const DMatrix& dv = qp_fcn_.output(it->qpf_def);
    copy(dv.begin(),dv.end(),it->step.begin());
  }
  
  // Expanded dual step (first part)
  if(!gauss_newton_){
    for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
      const DMatrix& dlam_v = qp_fcn_.output(it->qpf_defL);
      copy(dlam_v.begin(),dlam_v.end(),it->dlam.begin());
    }
  }
}

void SCPgenInternal::regularize(){
  casadi_assert(x_[0].n==2);
  
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
    std::cerr << "Regularization with " << reg_ << " to ensure positive definite Hessian." << endl;
    qpH_(0,0) += reg_;
    qpH_(1,1) += reg_;
  }
}

void SCPgenInternal::solve_qp(){  
  // Solve the QP
  qp_solver_.setInput(qpH_,QP_H);
  qp_solver_.setInput(qpG_,QP_G);
  qp_solver_.setInput(qpA_,QP_A);
  std::transform(lbu_.begin(),lbu_.end(), x_[0].opt.begin(),qp_solver_.input(QP_LBX).begin(),std::minus<double>());
  std::transform(ubu_.begin(),ubu_.end(), x_[0].opt.begin(),qp_solver_.input(QP_UBX).begin(),std::minus<double>()); 
  std::transform(lbg_.begin(),lbg_.end(), qpB_.begin(),qp_solver_.input(QP_LBA).begin(),std::minus<double>());
  std::transform(ubg_.begin(),ubg_.end(), qpB_.begin(),qp_solver_.input(QP_UBA).begin(),std::minus<double>());
  
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

void SCPgenInternal::line_search(int& ls_iter, bool& ls_success){
  // Make sure that we have a decent direction 
  if(!gauss_newton_){
    // Get the reduced Hessian times the step
    const vector<double>& qpH_times_du = exp_fcn_.output(exp_curve_).data();
    casadi_assert(qpH_times_du.size()==n_);
    
    // Scalar product with du to get gain
    double gain = 0;
    for(int i=0; i<n_; ++i){
      gain += x_[0].step[i] * qpH_times_du[i];
    }
    if (gain < 0){
      casadi_warning("Hessian indefinite in the search direction");
    }
  }

  // Calculate penalty parameter of merit function
  sigma_ = 0;
  sigma_ = std::max(sigma_,1.01*norm_inf(qp_solver_.output(QP_LAMBDA_X).data()));
  sigma_ = std::max(sigma_,1.01*norm_inf(qp_solver_.output(QP_LAMBDA_A).data()));
  
  // Calculate L1-merit function in the actual iterate
  double l1_infeas = primalInfeasibility();

  // Right-hand side of Armijo condition
  double F_sens = exp_fcn_.output(exp_osens_).toScalar();
  double L1dir = F_sens - sigma_ * l1_infeas;
  double L1merit = obj_k_ + sigma_ * l1_infeas;
  
  // Storing the actual merit function value in a list
  merit_mem_.push_back(L1merit);
  if (merit_mem_.size() > merit_memsize_){
    merit_mem_.pop_front();
  }

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
    for(vector<Var>::iterator it=x_.begin(); it!=x_.end(); ++it){
      for(int i=0; i<it->n; ++i){
	it->opt[i] += (t-t_prev) * it->step[i];
      }
    }
    
    // Take the dual step
    if(!gauss_newton_){
      for(int i=0; i<m_; ++i){
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

void SCPgenInternal::eval_exp(){
  // Pass current parameter guess
  exp_fcn_.setInput(x_[1].opt, exp_p_);

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

  // Expanded primal step (second part)
  for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
    const DMatrix& dv = exp_fcn_.output(it->exp_def);
    transform(dv.begin(),dv.end(),it->step.begin(),it->step.begin(),std::minus<double>());
  }
  
  // Expanded dual step (second part)
  if(!gauss_newton_){
    for(vector<Var>::iterator it=x_.begin()+1; it!=x_.end(); ++it){
      const DMatrix& dlam_v = exp_fcn_.output(it->exp_defL);
      transform(dlam_v.begin(),dlam_v.end(),it->dlam.begin(),it->dlam.begin(),std::minus<double>());
    }
  }  
}
  

} // namespace CasADi
