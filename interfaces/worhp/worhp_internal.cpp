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

#include "worhp_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include <ctime>

using namespace std;

namespace CasADi{

WorhpInternal::WorhpInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : NLPSolverInternal(F,G,H,J), GF_(GF){
  casadi_warning("This is just a boilerplate implementation - WorhpSolver does NOT work yet");
  

  
/**
  addOption("pass_nonlinear_variables", OT_BOOLEAN, true);
  addOption("print_time", OT_BOOLEAN, true, "print information about execution time");
  
  // Monitors
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_g|eval_jac_g|eval_grad_f", true);

  // Set pointers to zero
  app = 0;
  userclass = 0;

  // Start the application
  app = new Worhp::WorhpApplication();

  // Get all options available in WORHP
  map<string, Worhp::SmartPtr<Worhp::RegisteredOption> > regops = app->RegOptions()->RegisteredOptionsList();
  for(map<string, Worhp::SmartPtr<Worhp::RegisteredOption> >::const_iterator it=regops.begin(); it!=regops.end(); ++it){
    // Option identifier
    string opt_name = it->first;
    
    // Short description goes here, even though we do have a longer description
    string opt_desc = it->second->ShortDescription() + " (see WORHP documentation)";
    
    // Get the type
    Worhp::RegisteredOptionType worhp_type = it->second->Type();
    opt_type casadi_type;
    
    // Map Worhp option category to a CasADi options type
    switch(worhp_type){
      case Worhp::OT_Number:    casadi_type = OT_REAL;          break;
      case Worhp::OT_Integer:   casadi_type = OT_INTEGER;       break;
      case Worhp::OT_String:    casadi_type = OT_STRING;        break;
      case Worhp::OT_Unknown:   continue; // NOTE: No mechanism to handle OT_Unknown options
      default:                  continue; // NOTE: Unknown Worhp options category
    }
    
    addOption(opt_name, casadi_type, GenericType(), opt_desc);
    
    // Set default values of WORHP options 
    if (casadi_type == OT_REAL) {
      setDefault(opt_name,it->second->DefaultNumber());
    } else if (casadi_type == OT_INTEGER) {
      setDefault(opt_name,it->second->DefaultInteger());
    } else if (casadi_type == OT_STRING) {
      setDefault(opt_name,it->second->DefaultString());
    };
    
    // Save to map containing WORHP specific options
    ops_[opt_name] = casadi_type;
  }
  */
}


WorhpInternal::~WorhpInternal(){
/**
  if(app) delete app;

  // delete the smart pointer;
  if(userclass != 0){
    Worhp::SmartPtr<Worhp::TNLP> *ucptr = (Worhp::SmartPtr<Worhp::TNLP>*)userclass;
    delete ucptr;
  }
  */
}

void WorhpInternal::init(){

  NLPSolverInternal::init();
  
  worhp_o.initialised = false;
  worhp_w.initialised = false;
  worhp_p.initialised = false;
  worhp_c.initialised = false;
  
  worhp_o.n = n_;  // Number of variables
  worhp_o.m = m_;  // Number of constraints
  
  if (GF_.isNull()) GF_ = F_.jacobian();
  
  // Gradient of the objective function, remove?
  if(!GF_.isNull()) GF_.init();
  if(!GF_.isNull()) {
    casadi_assert_message(GF_.getNumInputs()>=1, "Wrong number of input arguments to GF");
    casadi_assert_message(GF_.getNumOutputs()>=1, "Wrong number of output arguments to GF");
    casadi_assert_message(GF_.input().numel()==n_,"Inconsistent dimensions");
    casadi_assert_message((GF_.output().size1()==n_ && GF_.output().size2()==1) || (GF_.output().size1()==1 && GF_.output().size2()==n_),"Inconsistent dimensions");
  }
  
  
  
  worhp_w.DF.nnz = GF_.output().size(); // Gradient of f
  worhp_w.DG.nnz = J_.output().size();  // Jacobian of G
  worhp_w.HM.nnz = worhp_o.n;  // Todo: what to do here?                 

  /* Data structure initialisation. */
  WorhpInit(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
  if (worhp_c.status != FirstCall) {
    casadi_error("Main: Initialisation failed.");
  }
  
  
  
}

void WorhpInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);

  // Set the static parameter
  if (!F_.isNull()) {
    if (F_.getNumInputs()==2) F_.setInput(input(NLP_P),1);
  }
  if (!G_.isNull()) {
    if (G_.getNumInputs()==2) G_.setInput(input(NLP_P),1);
  }
  if (!H_.isNull()) {
    if (H_.getNumInputs()==4) H_.setInput(input(NLP_P),1);
  }
  if (!J_.isNull()) {
    if (J_.getNumInputs()==2) J_.setInput(input(NLP_P),1);
  }
  if (!GF_.isNull()) {
    if (GF_.getNumInputs()==2) GF_.setInput(input(NLP_P),1);
  }
  
  input(NLP_X_INIT).getArray(worhp_o.X,n_);
  output(NLP_LAMBDA_X).getArray(worhp_o.Lambda,n_);
  input(NLP_LAMBDA_INIT).getArray(worhp_o.Mu,m_);
   
  input(NLP_LBX).getArray(worhp_o.XL,n_);
  input(NLP_UBX).getArray(worhp_o.XU,n_);
  input(NLP_LBG).getArray(worhp_o.GL,m_);
  input(NLP_UBG).getArray(worhp_o.GU,m_);


  if (worhp_w.DF.NeedStructure) {
    CRSSparsity GFT = trans(GF_.output()).sparsity();
    const std::vector<int> & col = GFT.col();
    for (int i=0;i<col.size();++i) worhp_w.DF.row[i] = col[i] + 1; // Index-1 based
  }
  
  if (worhp_w.DG.NeedStructure) {
    CRSSparsity JT = trans(J_.output()).sparsity();
    const std::vector<int> & col = JT.col();
    const std::vector<int> & row = JT.rowind();
    for (int i=0;i<col.size();++i) worhp_w.DG.row[i] = col[i] + 1;
    for (int i=0;i<row.size();++i) worhp_w.DG.col[i] = row[i] + 1;
  }
  
  std::cout << "clear" << std::endl;
  
  
  /**
  if (worhp_w.HM.NeedStructure) {
    int nz=0;
    vector<int> rowind,col;
    trans(H_.output()).sparsity().getSparsityCRS(rowind,col);
    for(int r=0; r<rowind.size()-1; ++r)
      for(int el=rowind[r]; el<rowind[r+1]; ++el){
       if(col[el]<=r){
          worhp_w.HM.col[nz] = r;
          worhp_w.HM.row[nz] = col[el];
          nz++;
       }
      }
  }
  */
  
 
  
  // Reverse Communication loop
  while(worhp_c.status < TerminateSuccess &&  worhp_c.status > TerminateError) {

    if (GetUserAction(&worhp_c, callWorhp)) {
      Worhp(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
    }

    if (GetUserAction(&worhp_c, iterOutput)) {
      IterationOutput(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
      DoneUserAction(&worhp_c, iterOutput);
    }

    if (GetUserAction(&worhp_c, evalF)) {
      eval_f(worhp_o.X, worhp_o.F);
      DoneUserAction(&worhp_c, evalF);
    }

    if (GetUserAction(&worhp_c, evalG)) {
      eval_g(worhp_o.X, worhp_o.G);
      DoneUserAction(&worhp_c, evalG);
    }

    if (GetUserAction(&worhp_c, evalDF)) {
      eval_grad_f(worhp_o.X, worhp_w.ScaleObj, worhp_w.DF.val);
      DoneUserAction(&worhp_c, evalDF);
    }

    if (GetUserAction(&worhp_c, evalDG)) {
      eval_jac_g(worhp_o.X,worhp_w.DG.val);
      DoneUserAction(&worhp_c, evalDG);
    }

    if (GetUserAction(&worhp_c, evalHM)) {
      eval_h(worhp_o.X, worhp_w.ScaleObj, worhp_o.Mu, worhp_w.HM.val);
      DoneUserAction(&worhp_c, evalHM);
    }
    
    if (GetUserAction(&worhp_c, fidif)) {
      WorhpFidif(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
    }

  }
  
  StatusMsg(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
  WorhpFree(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
 
}

bool WorhpInternal::intermediate_callback(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value, int iter, double inf_pr, double inf_du,double mu,double d_norm,double regularization_size,double alpha_du,double alpha_pr,int ls_trials) {
/*
  try {
    log("intermediate_callback started");
    double time1 = clock();
    if (!callback_.isNull()) {
#ifdef WITH_WORHP_CALLBACK 
      copy(x,x+n_,callback_.input(NLP_X_OPT).begin());
      
      vector<double>& lambda_x = callback_.input(NLP_LAMBDA_X).data();
      for(int i=0; i<lambda_x.size(); ++i){
        lambda_x[i] = z_U[i]-z_L[i];
      }
      copy(lambda,lambda+m_,callback_.input(NLP_LAMBDA_G).begin());
      copy(g,g+m_,callback_.input(NLP_G).begin());
      
#endif // WITH_WORHP_CALLBACK 

#ifndef WITH_WORHP_CALLBACK 
   if (iter==0) {
      cerr << "Warning: intermediate_callback is disfunctional in your installation. You will only be able to use getStats(). See https://sourceforge.net/apps/trac/casadi/wiki/enableWorhpCallback to enable it." << endl;
   }
#endif // WITH_WORHP_CALLBACK 
      callback_.input(NLP_COST).at(0) = obj_value;
      callback_->stats_["iter"] = iter;
      callback_->stats_["inf_pr"] = inf_pr;
      callback_->stats_["inf_du"] = inf_du;
      callback_->stats_["mu"] = mu;
      callback_->stats_["d_norm"] = d_norm;
      callback_->stats_["regularization_size"] = regularization_size;
      callback_->stats_["alpha_pr"] = alpha_pr;
      callback_->stats_["alpha_du"] = alpha_du;
      callback_->stats_["ls_trials"] = ls_trials;
      callback_.evaluate();
      double time2 = clock();
      t_callback_fun_ += double(time2-time1)/CLOCKS_PER_SEC;
      return  !callback_.output(0).at(0);
    } else {
      return 1;
    }
  } catch (exception& ex){
    if (getOption("iteration_callback_ignore_errors")) {
      cerr << "intermediate_callback: " << ex.what() << endl;
    } else {
      throw ex;
    }
  }
  */
}

void WorhpInternal::finalize_solution(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value){
/*
  try {
    // Get primal solution
    copy(x,x+n_,output(NLP_X_OPT).begin());

    // Get optimal cost
    output(NLP_COST).at(0) = obj_value;

    // Get dual solution (simple bounds)
    vector<double>& lambda_x = output(NLP_LAMBDA_X).data();
    for(int i=0; i<lambda_x.size(); ++i){
      lambda_x[i] = z_U[i]-z_L[i];
    }

    // Get dual solution (nonlinear bounds)
    copy(lambda,lambda+m_,output(NLP_LAMBDA_G).begin());
    
    // Get the constraints
    copy(g,g+m_,output(NLP_G).begin());
    
  } catch (exception& ex){
    cerr << "finalize_solution failed: " << ex.what() << endl;
  }
  */
}

bool WorhpInternal::eval_h(const double* x, double obj_factor, const double* lambda, double* values){
  try{
    log("eval_h started");
    double time1 = clock();
    // Number of inputs to the hessian
    int n_hess_in = H_.getNumInputs();
    
    // Pass input
    H_.setInput(x);
    if(n_hess_in>1){
      H_.setInput(lambda, n_hess_in==4? 2 : 1);
      H_.setInput(obj_factor, n_hess_in==4? 3 : 2);
    }

    // Evaluate
    H_.evaluate();

    // Scale objective
    if(n_hess_in==1 && obj_factor!=1.0){
      for(vector<double>::iterator it=H_.output().begin(); it!=H_.output().end(); ++it){
        *it *= obj_factor;
      }
    }

    // Get results
    H_.output().get(values,SPARSESYM);
      
    double time2 = clock();
    t_eval_h_ += double(time2-time1)/CLOCKS_PER_SEC;
    log("eval_h ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_h failed: " << ex.what() << endl;
    return false;
  }
}

bool WorhpInternal::eval_jac_g(const double* x,double* values){
  try{
    log("eval_jac_g started");
    
    // Quich finish if no constraints
    if(m_==0){
      log("eval_jac_g quick return (m==0)");
      return true;
    }
    
    double time1 = clock();

    // Pass the argument to the function
    J_.setInput(x);
    
     // Evaluate the function
    J_.evaluate();

    // Get the output
    J_.getOutput(values);
    
    if(monitored("eval_jac_g")){
      cout << "x = " << J_.input().data() << endl;
      cout << "J = " << endl;
      J_.output().printSparse();
    }
    
    double time2 = clock();
    t_eval_jac_g_ += double(time2-time1)/CLOCKS_PER_SEC;
    
    log("eval_jac_g ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_jac_g failed: " << ex.what() << endl;
    return false;
  }
}

bool WorhpInternal::eval_f(const double* x, double& obj_value)
{
  try {
    log("eval_f started");
    
    // Log time
    double time1 = clock();

    // Pass the argument to the function
    F_.setInput(x);
      
    // Evaluate the function
    F_.evaluate();

    // Get the result
    F_.getOutput(obj_value);

    // Printing
    if(monitored("eval_f")){
      cout << "x = " << F_.input() << endl;
      cout << "obj_value = " << obj_value << endl;
    }

    double time2 = clock();
    t_eval_f_ += double(time2-time1)/CLOCKS_PER_SEC;

    log("eval_f ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_f failed: " << ex.what() << endl;
    return false;
  }
}

bool WorhpInternal::eval_g(const double* x, double* g)
{
  try {
    log("eval_g started");
    double time1 = clock();

    // Pass the argument to the function
    G_.setInput(x);

    // Evaluate the function and tape
    G_.evaluate();

    // Ge the result
    G_.getOutput(g);

    // Printing
    if(monitored("eval_g")){
      cout << "x = " << G_.input() << endl;
      cout << "g = " << G_.output() << endl;
    }

      
    double time2 = clock();
    t_eval_g_ += double(time2-time1)/CLOCKS_PER_SEC;
    
    log("eval_g ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_g failed: " << ex.what() << endl;
    return false;
  }
}

bool WorhpInternal::eval_grad_f(const double* x,double scale , double* grad_f )
{
  try {
    log("eval_grad_f started");
    double time1 = clock();
    
    // If no gradient function has been provided, use AD adjoint
    if(GF_.isNull()){
    
      // Pass the argument to the function
      F_.setInput(x);
      
      // Give a seed to the function
      F_.setAdjSeed(scale);

      // Evaluate, adjoint mode
      F_.evaluate(0,1);

      // Get the result
      F_.getAdjSens(grad_f);

      // Printing
      if(monitored("eval_grad_f")){
        cout << "grad_f = " << F_.adjSens() << endl;
      }
      
    } else {
      
      // Pass the argument to the function
      GF_.setInput(x);
      
      // Evaluate, adjoint mode
      GF_.evaluate();

      GF_.output()*=scale;
      // Get the result
      GF_.getOutput(grad_f);
      
      // Printing
      if(monitored("eval_grad_f")){
        cout << "grad_f = " << GF_.output() << endl;
      }
    }
    
    double time2 = clock();
    t_eval_grad_f_ += double(time2-time1)/CLOCKS_PER_SEC;

    // Check the result for regularity
    for(int i=0; i<n_; ++i){
        if(isnan(grad_f[i]) || isinf(grad_f[i])){
          log("eval_grad_f: result not regular");
          return false;
      }
    }

    log("eval_grad_f ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_jac_f failed: " << ex.what() << endl;
    return false;
  }
}


void WorhpInternal::get_nlp_info(int& n, int& m, int& nnz_jac_g,int& nnz_h_lag)
{
/*
  try {
    n = n_;               // number of variables
    m = m_;               // number of constraints

    // Get Jacobian sparsity pattern
    if(G_.isNull())
      nnz_jac_g = 0;
    else
      nnz_jac_g = J_.output().size();

    // Get Hessian sparsity pattern
    if(exact_hessian_)
      nnz_h_lag = H_.output().sparsity().sizeL();
    else
      nnz_h_lag = 0;
  } catch (exception& ex){
    cerr << "get_nlp_info failed: " << ex.what() << endl;
  }
  */
}

int WorhpInternal::get_number_of_nonlinear_variables() const{
/*
  try {
    if(H_.isNull() || !bool(getOption("pass_nonlinear_variables"))){
      // No Hessian has been interfaced
      return -1;
    } else {
      // Number of variables that appear nonlinearily
      int nv = 0;
      
      // Loop over the rows
      for(int i=0; i<H_.output().size1(); ++i){
        // If the row contains any non-zeros, the corresponding variable appears nonlinearily
        if(H_.output().rowind(i)!=H_.output().rowind(i+1))
          nv++;
      }
      
      // Return the number
      return nv;
    }
  } catch (exception& ex){
    cerr << "get_number_of_nonlinear_variables failed: " << ex.what() << endl;
    return -1;
  }
  */
}

bool WorhpInternal::get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars) const{
/*
  try {
    // Running index
    int el = 0;
    
    // Loop over the rows
    for(int i=0; i<H_.output().size1(); ++i){
      // If the row contains any non-zeros, the corresponding variable appears nonlinearily
      if(H_.output().rowind(i)!=H_.output().rowind(i+1)){
        pos_nonlin_vars[el++] = i;
      }
    }
    
    // Assert number and return
    casadi_assert(el==num_nonlin_vars);
    return true;
  } catch (exception& ex){
    cerr << "get_list_of_nonlinear_variables failed: " << ex.what() << endl;
    return false;
  }
  */
}

} // namespace CasADi
