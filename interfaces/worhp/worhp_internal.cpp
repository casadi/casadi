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
#include "casadi/mx/mx_tools.hpp"
#include "casadi/matrix/sparsity_tools.hpp"
#include "casadi/fx/mx_function.hpp"
#include <ctime>

using namespace std;

namespace CasADi{

WorhpInternal::WorhpInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : NLPSolverInternal(F,G,H,J), GF_(GF){
  casadi_warning("WorhpSolver is not mature yet");

}


WorhpInternal::~WorhpInternal(){

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
  
  
  // Worhp uses the CS format internally, hence it is the preferred sparse matrix format.
  
  worhp_w.DF.nnz = GF_.output().size(); // Gradient of f
  worhp_w.DG.nnz = J_.output().size();  // Jacobian of G
  
  
  if (hasSetOption("generate_hessian") && getOption("generate_hessian")) {
    std::vector< MX > input = H_.symbolicInput();
    MX H = H_.call(input)[0];
    H = vertcat(vec(H(lowerSparsity(H.sparsity(),false))),vec(H(sp_diag(n_))));
    
    H_tril_ = MXFunction(input,H);
    H_tril_.init();
    
    worhp_w.HM.nnz = H_tril_.output().size();
  }

         
  /* Data structure initialisation. */
  WorhpInit(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
  if (worhp_c.status != FirstCall) {
    casadi_error("Main: Initialisation failed.");
  }
  
  if (worhp_w.DF.NeedStructure) {
    vector<int> row,col;
    GF_.output().sparsity().getSparsity(col,row); // transpose
    for (int i=0;i<row.size();++i) worhp_w.DF.row[i] = row[i] + 1; // Index-1 based
  }
  
  if (worhp_w.DG.NeedStructure) {

    
    vector<int> row,col;
    trans(J_.output()).sparsity().getSparsity(col,row);
    std::vector< MX > J = J_.symbolicInput();
      
    casadi_assert(col.size()==worhp_w.DG.nnz);
    casadi_assert(row.size()==worhp_w.DG.nnz);
    
    for (int i=0;i<col.size();++i) worhp_w.DG.col[i] = col[i] + 1;
    for (int i=0;i<row.size();++i) worhp_w.DG.row[i] = row[i] + 1;
    

  }
  
    

  if (hasSetOption("generate_hessian") && getOption("generate_hessian")) {
    log("generate_hessian sparsity");
    int nz=0;
    if (worhp_w.HM.NeedStructure) {
      vector<int> row,col;
      lowerSparsity(H_.output().sparsity(),false).getSparsity(row,col);

      for (int i=0;i<col.size();++i) worhp_w.HM.col[i] = col[i] + 1;
      for (int i=0;i<row.size();++i) worhp_w.HM.row[i] = row[i] + 1;
      
      vector<int> rowd,cold;
      H_.output()(sp_diag(n_)).sparsity().getSparsity(rowd,cold);
      
      casadi_assert(worhp_w.HM.nnz<=worhp_w.HM.dim_row);
      casadi_assert(worhp_w.HM.nnz<=worhp_w.HM.dim_col);
      
      casadi_assert(cold.size()+col.size()==worhp_w.HM.nnz);
      casadi_assert(rowd.size()+row.size()==worhp_w.HM.nnz);
      
      for (int i=0;i<cold.size();++i) worhp_w.HM.col[i+col.size()] = cold[i] + 1;
      for (int i=0;i<rowd.size();++i) worhp_w.HM.row[i+row.size()] = rowd[i] + 1;

    }
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
  if (!H_tril_.isNull()) {
    if (H_tril_.getNumInputs()==4) H_tril_.setInput(input(NLP_P),1);
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

  
  // Reverse Communication loop
  while(worhp_c.status < TerminateSuccess &&  worhp_c.status > TerminateError) {

    if (GetUserAction(&worhp_c, callWorhp)) {
      Worhp(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
    }

    if (GetUserAction(&worhp_c, iterOutput)) {
      if (!callback_.isNull()) {
        // Copy outputs
        copy(worhp_o.X,worhp_o.X+n_,callback_.input(NLP_X_OPT).begin());
        callback_.input(NLP_COST)[0] = worhp_o.F;
        copy(worhp_o.G,worhp_o.G+n_,callback_.input(NLP_G).begin());
        copy(worhp_o.Lambda,worhp_o.Lambda+n_,callback_.input(NLP_LAMBDA_X).begin());
        copy(worhp_o.Mu,worhp_o.Mu+n_,callback_.input(NLP_LAMBDA_G).begin());
        
        callback_.evaluate();
      }
    
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
  
  // Copy outputs
  copy(worhp_o.X,worhp_o.X+n_,output(NLP_X_OPT).begin());
  output(NLP_COST)[0] = worhp_o.F;
  copy(worhp_o.G,worhp_o.G+n_,output(NLP_G).begin());
  copy(worhp_o.Lambda,worhp_o.Lambda+n_,output(NLP_LAMBDA_X).begin());
  copy(worhp_o.Mu,worhp_o.Mu+n_,output(NLP_LAMBDA_G).begin());
  
  StatusMsg(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
  WorhpFree(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
 
}

bool WorhpInternal::eval_h(const double* x, double obj_factor, const double* lambda, double* values){
  try{
    log("eval_h started");
    double time1 = clock();
    // Number of inputs to the hessian
    int n_hess_in = H_tril_.getNumInputs();
    
    // Pass input
    H_tril_.setInput(x);
    if(n_hess_in>1){
      H_tril_.setInput(lambda, n_hess_in==4? 2 : 1);
      H_tril_.setInput(obj_factor, n_hess_in==4? 3 : 2);
    }

    // Evaluate
    H_tril_.evaluate();

    // Scale objective
    if(n_hess_in==1 && obj_factor!=1.0){
      for(vector<double>::iterator it=H_tril_.output().begin(); it!=H_tril_.output().end(); ++it){
        *it *= obj_factor;
      }
    }

    // Get results
    H_tril_.output().get(values);
      
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
    trans(J_.output()).get(values);
    
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


} // namespace CasADi
