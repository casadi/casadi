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

  SCPgenInternal::SCPgenInternal(const FX& nlp) : NLPSolverInternal(nlp){
    casadi_warning("SCPgen is under development");
    addOption("qp_solver",         OT_QPSOLVER,   GenericType(),    "The QP solver to be used by the SQP method");
    addOption("qp_solver_options", OT_DICTIONARY, GenericType(),    "Options to be passed to the QP solver");
    addOption("hessian_approximation", OT_STRING, "exact",          "gauss-newton|exact");
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
    string compiler = getOption("compiler");
    gauss_newton_ = getOption("hessian_approximation") == "gauss-newton";
    if(gauss_newton_){
      casadi_assert(nlp_.output(NLP_F).size()>1);
    } else {
      casadi_assert(nlp_.output(NLP_F).size()==1);
    }

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

    MXFunction fg = shared_cast<MXFunction>(nlp_);
    if(fg.isNull()){
      vector<MX> nlp_in = nlp_.symbolicInput();
      vector<MX> nlp_out = nlp_.call(nlp_in);
      fg = MXFunction(nlp_in,nlp_out);
    }
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
    g_lam_.resize(ng_,0);
    g_dlam_.resize(ng_,0);
    qpH_times_du_.resize(nx_);
 
    x_init_.resize(nx_,0);
    x_opt_.resize(nx_,0);
    x_step_.resize(nx_,0);
    x_lam_.resize(nx_,0);
    x_dlam_.resize(nx_,0);

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
    MX f;

    // Multipliers
    MX g_lam;

    // Definition of the lifted dual variables
    MX p_defL, gL_defL;

    if(gauss_newton_){    
      // Least square objective
      f = inner_prod(vdef_out[0],vdef_out[0])/2;
      gL_defL = vdef_out[0];
      b_gn_.resize(gL_defL.size(),numeric_limits<double>::quiet_NaN());

    } else {
      // Scalar objective function
      f = vdef_out[0];
    
      // Lagrange multipliers corresponding to the definition of the dependent variables
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
      vdef_fcn.eval(vdef_in,vdef_out,fseed,fsens,aseed,asens);
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

      if(verbose_){
        cout << "Generated the gradient of the Lagrangian." << endl;
      }
    }
    gL_.resize(nx_, numeric_limits<double>::quiet_NaN());
    gf_.resize(nx_, numeric_limits<double>::quiet_NaN());

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
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        res_fcn_in.push_back(it->v_lam);  it->res_lam = n++;
      }
    }

    // Outputs
    vector<MX> res_fcn_out;
    n=0;
    res_fcn_out.push_back(f);                              res_f_ = n++;
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
      res_fcn_ = dynamicCompilation(res_fcn,"res_fcn","residual function",compiler);
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
    ex[0] = f;
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

    MX f_z = ex[0];
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
  
    // Add multipliers to function inputs
    if(!gauss_newton_){
      n = mfcn_in.size();
      mfcn_in.push_back(g_lam);                          mod_g_lam_ = n++;
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        mfcn_in.push_back(it->d_lam);                    it->mod_lam = n++;
      }
    }

    // Add gradient of the Lagrangian
    n = mfcn_out.size();
    mfcn_out.push_back(f_z);                             mod_f_ = n++;
    mfcn_out.push_back(gL_z);                            mod_gl_ = n++;  
  
    // Lagrangian gradient function
    MXFunction lgrad(mfcn_in,mfcn_out);
    lgrad.init();
  
    // Jacobian of the constraints
    MX jac = lgrad.jac(mod_x_,mod_g_);
    log("Formed Jacobian of the constraints.");

    // Hessian of the Lagrangian
    MX hes = lgrad.jac(mod_x_,mod_gl_,false,!gauss_newton_);
    if(gauss_newton_){
      log("Formed square root of Gauss-Newton Hessian.");
    } else {
      log("Formed Hessian of the Lagrangian.");
    }

    // Matrices in the reduced QP
    n=0;
    vector<MX> mat_out;  
    mat_out.push_back(jac);                             mat_jac_ = n++;
    mat_out.push_back(hes);                             mat_hes_ = n++;
    MXFunction mat_fcn(mfcn_in,mat_out);
    mat_fcn.init();
  
    // Generate c code and load as DLL
    if(codegen_){
      mat_fcn_ = dynamicCompilation(mat_fcn,"mat_fcn","Matrices function",compiler);
    } else {
      mat_fcn_ = mat_fcn;
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
    mfcn.eval(mfcn_in,mfcn_out,mfcn_fwdSeed,mfcn_fwdSens,mfcn_adjSeed,mfcn_adjSens);
  
    // Vector(s) b in Lifted Newton
    MX b_gf = mfcn_fwdSens[0][mod_gl_];
    MX b_g = mfcn_fwdSens[0][mod_g_];
  
    // Make sure that the vectors are dense
    makeDense(b_gf);
    makeDense(b_g);
  
    // Tangent function
    vector<MX> vec_fcn_out;
    n=0;
    vec_fcn_out.push_back(b_gf);                              vec_gf_ = n++;
    vec_fcn_out.push_back(b_g);                               vec_g_ = n++;  
    casadi_assert(n==vec_fcn_out.size());
  
    MXFunction vec_fcn(mfcn_in,vec_fcn_out);
    vec_fcn.setOption("name","vec_fcn");
    vec_fcn.setOption("number_of_fwd_dir",0);
    vec_fcn.setOption("number_of_adj_dir",0);
    vec_fcn.init();
    if(verbose_){
      cout << "Generated linearization function ( " << vec_fcn.getAlgorithmSize() << " nodes)." << endl;
    }
  
    // Generate c code and load as DLL
    if(codegen_){
      vec_fcn_ = dynamicCompilation(vec_fcn,"vec_fcn","linearization function",compiler);
    } else {
      vec_fcn_ = vec_fcn;
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
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      mfcn_fwdSeed[0][it->mod_var] = -it->d;
    }
    if(!gauss_newton_){
      mfcn_fwdSeed[0][mod_g_lam_] = g_dlam;
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        mfcn_fwdSeed[0][it->mod_lam] = -it->d_lam;
      }
    }

    mfcn.eval(mfcn_in,mfcn_out,mfcn_fwdSeed,mfcn_fwdSens,mfcn_adjSeed,mfcn_adjSens);    
  
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
      exp_fcn_ = dynamicCompilation(exp_fcn,"exp_fcn","step expansion function",compiler);
    } else {
      exp_fcn_ = exp_fcn;
    }  
  
    // Allocate QP data
    CRSSparsity sp_tr_B_obj = mat_fcn_.output(mat_hes_).sparsity().transpose();
    qpH_ = DMatrix(sp_tr_B_obj.patternProduct(sp_tr_B_obj));
    qpA_ = mat_fcn_.output(mat_jac_);
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
    const vector<double>& x_init = input(NLP_SOLVER_X0).data();
    const vector<double>& lbx = input(NLP_SOLVER_LBX).data();
    const vector<double>& ubx = input(NLP_SOLVER_UBX).data();
    const vector<double>& lbg = input(NLP_SOLVER_LBG).data();
    const vector<double>& ubg = input(NLP_SOLVER_UBG).data();  
  
    copy(x_init.begin(),x_init.end(),x_init_.begin());
    copy(lbx.begin(),lbx.end(),x_lb_.begin());
    copy(ubx.begin(),ubx.end(),x_ub_.begin());
    copy(lbg.begin(),lbg.end(),g_lb_.begin());
    copy(ubg.begin(),ubg.end(),g_ub_.begin());
  
    if(v_.size()>0){
      // Initialize lifted variables using the generated function
      vinit_fcn_.setInput(x_init_,0);
      vinit_fcn_.setInput(input(NLP_SOLVER_P),1);
      vinit_fcn_.evaluate();    
      for(int i=0; i<v_.size(); ++i){
        vinit_fcn_.getOutput(v_[i].init,i);
      }
    }
    if(verbose_){
      cout << "Passed initial guess" << endl;
    }

    // Reset dual guess
    fill(g_lam_.begin(),g_lam_.end(),0);
    fill(g_dlam_.begin(),g_dlam_.end(),0);    
    fill(x_lam_.begin(),x_lam_.end(),0);
    fill(x_dlam_.begin(),x_dlam_.end(),0);
    if(!gauss_newton_){
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        fill(it->lam.begin(),it->lam.end(),0);
        fill(it->dlam.begin(),it->dlam.end(),0);
      }
    }
  
    // Objective value
    f_ = numeric_limits<double>::quiet_NaN();

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
    t_eval_mat_ = t_eval_res_ = t_eval_vec_ = t_eval_exp_ = t_solve_qp_ = 0;

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

      // Evaluate the vectors in the condensed QP
      eval_vec();

      // Evaluate the matrices in the condensed QP
      eval_mat();

      // 1-norm of the primal infeasibility
      double pr_inf = primalInfeasibility();
    
      // 1-norm of the dual infeasibility
      double du_inf = dualInfeasibility();
    
      // Print header occasionally
      if(iter % 10 == 0) printIteration(cout);
    
      // Printing information about the actual iterate
      printIteration(cout,iter,f_,pr_inf,du_inf,reg_,ls_iter,ls_success);

      // Checking convergence criteria
      bool converged = pr_inf <= tol_pr_ && pr_step_ <= tol_pr_step_ && reg_ <= tol_reg_;
      converged = converged && du_inf <= tol_du_;
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
      if(f_!=f_ || pr_step_ != pr_step_ || pr_inf != pr_inf){
        cout << "CasADi::SCPgen: Aborted, nan detected" << endl;
        break;
      }
    
      // Start a new iteration
      iter++;
        
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
    cout << "optimal cost = " << f_ << endl;

    // Save results to outputs
    output(NLP_SOLVER_F).set(f_);
    output(NLP_SOLVER_X).set(x_opt_);
    output(NLP_SOLVER_LAM_G).set(g_lam_);
    output(NLP_SOLVER_LAM_X).set(x_lam_);
    output(NLP_SOLVER_G).set(g_);
  
    // Write timers
    if(print_time_){
      cout << endl;
      cout << "time spent in eval_mat:    " << setw(9) << t_eval_mat_ << " s." << endl;
      cout << "time spent in eval_res:    " << setw(9) << t_eval_res_ << " s." << endl;
      cout << "time spent in eval_vec:    " << setw(9) << t_eval_vec_ << " s." << endl;
      cout << "time spent in eval_exp:    " << setw(9) << t_eval_exp_ << " s." << endl;
      cout << "time spent in solve_qp:    " << setw(9) << t_solve_qp_ << " s." << endl;
      cout << "time spent in main loop:   " << setw(9) << t_mainloop_ << " s." << endl;
    }

    // Save statistics
    stats_["iter_count"] = iter;

    cout << endl;
  }  

  double SCPgenInternal::primalInfeasibility(){
    // L1-norm of the primal infeasibility
    double pr_inf = 0;
  
    // Simple bounds
    for(int i=0; i<nx_; ++i) pr_inf +=  std::max(x_opt_[i]-x_ub_[i],0.);
    for(int i=0; i<nx_; ++i) pr_inf +=  std::max(x_lb_[i]-x_opt_[i],0.);

    // Lifted variables
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      for(int i=0; i<it->n; ++i) pr_inf += ::fabs(it->res[i]);
    }
  
    // Nonlinear bounds
    for(int i=0; i<ng_; ++i) pr_inf += std::max(g_[i]-g_ub_[i],0.);
    for(int i=0; i<ng_; ++i) pr_inf += std::max(g_lb_[i]-g_[i],0.);
  
    return pr_inf;
  }  

  double SCPgenInternal::dualInfeasibility(){

    // L1-norm of the dual infeasibility
    double du_inf = 0;
  
    // Lifted variables
    for(int i=0; i<nx_; ++i) du_inf += ::fabs(gL_[i]);

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
    stream << setprecision(2) << du_inf;
    stream << setw(11) << setprecision(2) << pr_step_;
    stream << setw(11);
    stream << setprecision(2) << du_step_;
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

  void SCPgenInternal::eval_mat(){
    // Get current time
    double time1 = clock();

    // Pass parameters
    mat_fcn_.setInput(input(NLP_SOLVER_P),mod_p_);

    // Pass primal step/variables  
    mat_fcn_.setInput(x_opt_, mod_x_);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      mat_fcn_.setInput(it->res, it->mod_var);
    }

    // Pass dual steps/variables
    if(!gauss_newton_){
      mat_fcn_.setInput(g_lam_,mod_g_lam_);
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        mat_fcn_.setInput(it->resL, it->mod_lam);
      }
    }

    // Evaluate condensed Hessian
    mat_fcn_.evaluate();

    // Get the Jacobian
    mat_fcn_.getOutput(qpA_,mat_jac_);

    if(gauss_newton_){
      // Gauss-Newton Hessian
      const DMatrix& B_obj =  mat_fcn_.output(mat_hes_);
      fill(qpH_.begin(),qpH_.end(),0);
      DMatrix::mul_no_alloc_tn(B_obj,B_obj,qpH_);

      // Gradient of the objective in Gauss-Newton
      fill(gf_.begin(),gf_.end(),0);
      DMatrix::mul_no_alloc_tn(B_obj,b_gn_,gf_);
    } else {
      // Exact Hessian
      mat_fcn_.getOutput(qpH_,mat_hes_);
    }

    // Calculate the gradient of the lagrangian
    const vector<double> &qpA_data = qpA_.data();
    const vector<int> &qpA_rowind = qpA_.rowind();
    const vector<int> &qpA_col = qpA_.col();
    for(int i=0; i<nx_; ++i)  gL_[i] = gf_[i] + x_lam_[i];
    for(int i=0; i<ng_; ++i){
      for(int el=qpA_rowind[i]; el<qpA_rowind[i+1]; ++el){
        gL_[qpA_col[el]] += qpA_data[el]*g_lam_[i];
      }
    }

    double time2 = clock();
    t_eval_mat_ += double(time2-time1)/CLOCKS_PER_SEC;
  }

  void SCPgenInternal::eval_res(){
    // Get current time
    double time1 = clock();

    // Pass parameters
    res_fcn_.setInput(input(NLP_SOLVER_P),res_p_);

    // Pass primal variables to the residual function for initial evaluation
    res_fcn_.setInput(x_opt_,res_x_);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      res_fcn_.setInput(it->opt,it->res_var);
    }

    // Pass dual variables to the residual function for initial evaluation
    if(!gauss_newton_){
      res_fcn_.setInput(0.0,res_g_lam_);
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        res_fcn_.setInput(it->lam,it->res_lam);
      }
    }
  
    // Evaluate residual function
    res_fcn_.evaluate();

    // Get objective
    f_ = res_fcn_.output(res_f_).toScalar();

    // Get objective gradient
    if(gauss_newton_){
      res_fcn_.getOutput(b_gn_,res_gl_);
    } else {
      res_fcn_.getOutput(gf_,res_gl_);
    }

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
    res_fcn_.getOutput(output(NLP_SOLVER_LAM_P),res_p_d_);

    double time2 = clock();
    t_eval_res_ += double(time2-time1)/CLOCKS_PER_SEC;
  }

  void SCPgenInternal::eval_vec(){
    // Get current time
    double time1 = clock();

    // Pass current parameter guess
    vec_fcn_.setInput(input(NLP_SOLVER_P),mod_p_);

    // Pass primal step/variables
    vec_fcn_.setInput(x_opt_, mod_x_);
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      vec_fcn_.setInput(it->res, it->mod_var);
    }

    // Pass dual steps/variables
    if(!gauss_newton_){
      vec_fcn_.setInput(0.0,mod_g_lam_);
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        vec_fcn_.setInput(it->resL, it->mod_lam);
      }
    }

    // Evaluate to get QP
    vec_fcn_.evaluate();

    // Linear offset in the reduced QP
    transform(g_.begin(),g_.end(),vec_fcn_.output(vec_g_).begin(),qpB_.begin(),std::minus<double>());

    // Gradient of the objective in the reduced QP
    if(gauss_newton_){
      transform(b_gn_.begin(),b_gn_.end(),vec_fcn_.output(vec_gf_).begin(),b_gn_.begin(),std::minus<double>());
    } else {
      transform(gf_.begin(),gf_.end(),vec_fcn_.output(vec_gf_).begin(),gf_.begin(),std::minus<double>());
    }
  
    double time2 = clock();
    t_eval_vec_ += double(time2-time1)/CLOCKS_PER_SEC;
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
    qp_solver_.setInput(gf_,QP_G);
    qp_solver_.setInput(qpA_,QP_A);
    std::transform(x_lb_.begin(),x_lb_.end(), x_opt_.begin(),qp_solver_.input(QP_LBX).begin(),std::minus<double>());
    std::transform(x_ub_.begin(),x_ub_.end(), x_opt_.begin(),qp_solver_.input(QP_UBX).begin(),std::minus<double>()); 
    std::transform(g_lb_.begin(),g_lb_.end(), qpB_.begin(),qp_solver_.input(QP_LBA).begin(),std::minus<double>());
    std::transform(g_ub_.begin(),g_ub_.end(), qpB_.begin(),qp_solver_.input(QP_UBA).begin(),std::minus<double>());
  
    qp_solver_.evaluate();
  
    // Condensed primal step
    const DMatrix& du = qp_solver_.output(QP_PRIMAL);
    copy(du.begin(),du.end(),x_step_.begin());
  
    // Condensed dual step (simple bounds)
    const DMatrix& lam_x_new = qp_solver_.output(QP_LAMBDA_X);
    copy(lam_x_new.begin(),lam_x_new.end(),x_dlam_.begin());
    std::transform(x_dlam_.begin(),x_dlam_.end(),x_lam_.begin(),x_dlam_.begin(),std::minus<double>());

    // Condensed dual step (nonlinear bounds)
    const DMatrix& lam_g_new = qp_solver_.output(QP_LAMBDA_A);
    copy(lam_g_new.begin(),lam_g_new.end(),g_dlam_.begin());
    std::transform(g_dlam_.begin(),g_dlam_.end(),g_lam_.begin(),g_dlam_.begin(),std::minus<double>());

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
    for(int i=0; i<nx_; ++i) F_sens += x_step_[i] * gf_[i];
    double L1dir = F_sens - sigma_ * l1_infeas;
    double L1merit = f_ + sigma_ * l1_infeas;
  
    // Storing the actual merit function value in a list
    merit_mem_[merit_ind_] = L1merit;
    ++merit_ind_ %= merit_memsize_;

    // Stepsize
    double t = 1.0, t_prev = 0.0;
  
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
      for(int i=0; i<ng_; ++i)         g_lam_[i] += (t-t_prev) * g_dlam_[i];
      for(int i=0; i<nx_; ++i)  x_lam_[i] += (t-t_prev) * x_dlam_[i];
      if(!gauss_newton_){
        for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
          for(int i=0; i<it->n; ++i) it->lam[i] += (t-t_prev) * it->dlam[i];
        }
      }
    
      // Evaluate residual function to get objective and constraints (and residuals for the next iteration)
      eval_res();
      ls_iter++;      
    
      // Calculating merit-function in candidate
      l1_infeas = primalInfeasibility();
      L1merit_cand = f_ + sigma_ * l1_infeas;
    
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

  void SCPgenInternal::eval_exp(){
    // Get current time
    double time1 = clock();

    // Pass current parameter guess
    exp_fcn_.setInput(input(NLP_SOLVER_P), mod_p_);

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
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        exp_fcn_.setInput(it->resL,it->mod_lam);
      }
    }

    // Perform the step expansion
    exp_fcn_.evaluate();

    // Expanded primal step
    for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
      const DMatrix& dv = exp_fcn_.output(it->exp_def);
      copy(dv.begin(),dv.end(),it->step.begin());
    }
  
    // Expanded dual step
    if(!gauss_newton_){
      for(vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it){
        const DMatrix& dlam_v = exp_fcn_.output(it->exp_defL);
        copy(dlam_v.begin(),dlam_v.end(),it->dlam.begin());
      }
    }

    double time2 = clock();
    t_eval_exp_ += double(time2-time1)/CLOCKS_PER_SEC;
  }
  

} // namespace CasADi
