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

#include "knitro_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>

using namespace std;
namespace CasADi{

KnitroInternal::KnitroInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : F_(F), G_(G), H_(H), J_(J), GF_(GF){
  casadi_warning("KnitroInternal: the KNITRO interface is still experimental, more tests are needed");
  kc_handle_ = 0;
  
  addOption("contype", OT_INTEGERVECTOR);
}


KnitroInternal::~KnitroInternal(){
  // Free KNITRO memory
  if(kc_handle_){
/*    KTR_free(&kc_handle_);
    kc_handle_ = 0;*/
  }
}

void KnitroInternal::init(){
  // Initialize the functions
  F_.init();
  if(!G_.isNull()) G_.init();
  n_ = F_.input(0).numel();
  m_ = G_.isNull() ? 0 : G_.output(0).numel();

  // Call the init method of the base class
  NLPSolverInternal::init();

  // Create a Jacobian if it does not already exists
  if(!G_.isNull() && J_.isNull()){
    J_ = G_.jacobian();
  }
  if(!J_.isNull()) J_.init();
  if(!H_.isNull()) H_.init();
  if(!GF_.isNull()) GF_.init();

  // Commented out since I have not found out how to change the bounds
  // Allocate KNITRO memory block
/*  casadi_assert(kc_handle_==0);
  kc_handle_ = KTR_new();*/
  
}

void KnitroInternal::evaluate(int nfdir, int nadir){
  // Allocate KNITRO memory block (move back to init!)
  casadi_assert(kc_handle_==0);
  kc_handle_ = KTR_new();
  
  casadi_assert(nfdir==0 && nadir==0);
  casadi_assert(kc_handle_!=0);
  int status;
  
  // Jacobian sparsity
  vector<int> Jcol = J_.output().col();
  vector<int> Jrow = J_.output().sparsity().getRow();

  // Hessian sparsity
  int nnzH = H_.isNull() ? 0 : H_.output().sizeL();
  vector<int> Hcol(nnzH), Hrow(nnzH);
  if(nnzH>0){
    const vector<int> &rowind = H_.output().rowind();
    const vector<int> &col = H_.output().col();
    int nz=0;
    for(int r=0; r<rowind.size()-1; ++r){
      for(int el=rowind[r]; el<rowind[r+1] && col[el]<=r; ++el){
        Hcol[nz] = r;
        Hrow[nz] = col[el];
        nz++;
      }
    }
    casadi_assert(nz==nnzH);
    
    status = KTR_set_int_param_by_name(kc_handle_, "hessopt", KTR_HESSOPT_EXACT);
    casadi_assert_message(status==0, "KTR_set_int_param failed");
  } else {
    status = KTR_set_int_param_by_name(kc_handle_, "hessopt", KTR_HESSOPT_LBFGS);
    casadi_assert_message(status==0, "KTR_set_int_param failed");
  }

  // Set user set options
  for(std::map<std::string, double>::iterator it=double_param_.begin(); it!=double_param_.end(); ++it){
    status = KTR_set_double_param_by_name(kc_handle_, it->first.c_str(), it->second);
    if(status!=0){
      throw CasadiException("KnitroInternal::evaluate: cannot set " + it->first);
    }
  }
  
  for(std::map<std::string, int>::iterator it=int_param_.begin(); it!=int_param_.end(); ++it){
    status = KTR_set_int_param_by_name(kc_handle_, it->first.c_str(), it->second);
    if(status!=0){
      throw CasadiException("KnitroInternal::evaluate: cannot set " + it->first);
    }
  }

  // Type of constraints
  vector<int> cType(m_,KTR_CONTYPE_GENERAL);
  if(hasSetOption("contype")){
    vector<int> contype = getOption("contype");
    casadi_assert(contype.size()==cType.size());
    copy(contype.begin(),contype.end(),cType.begin());
  }
  
  // "Correct" upper and lower bounds
  for(vector<double>::iterator it=input(NLP_LBX).begin(); it!=input(NLP_LBX).end(); ++it)
    if(isinf(*it)) *it = -KTR_INFBOUND;
  for(vector<double>::iterator it=input(NLP_UBX).begin(); it!=input(NLP_UBX).end(); ++it)
    if(isinf(*it)) *it =  KTR_INFBOUND;
  for(vector<double>::iterator it=input(NLP_LBG).begin(); it!=input(NLP_LBG).end(); ++it)
    if(isinf(*it)) *it = -KTR_INFBOUND;
  for(vector<double>::iterator it=input(NLP_UBG).begin(); it!=input(NLP_UBG).end(); ++it)
    if(isinf(*it)) *it =  KTR_INFBOUND;
  
  // Initialize KNITRO
  status = KTR_init_problem(kc_handle_, n_, KTR_OBJGOAL_MINIMIZE, KTR_OBJTYPE_GENERAL,
                              &input(NLP_LBX).front(), &input(NLP_UBX).front(),
                              m_, &cType.front(), &input(NLP_LBG).front(), &input(NLP_UBG).front(),
                              Jcol.size(), &Jcol.front(), &Jrow.front(),
                              nnzH,
                              nnzH==0 ? 0 : &Hrow[0],
                              nnzH==0 ? 0 : &Hcol[0],
                              &input(NLP_X_INIT).front(),
                              0); // initial lambda
  casadi_assert_message(status==0, "KTR_init_problem failed");
  
  // Register callback functions
  status = KTR_set_func_callback(kc_handle_, &callback);
  casadi_assert_message(status==0, "KTR_set_func_callback failed");
  
  status = KTR_set_grad_callback(kc_handle_, &callback);
  casadi_assert_message(status==0, "KTR_set_grad_callbackfailed");
  
  if(nnzH>0){
    status = KTR_set_hess_callback(kc_handle_, &callback);
    casadi_assert_message(status==0, "KTR_set_hess_callbackfailed");
  }

  // Lagrange multipliers
  vector<double> lambda(n_+m_);

  // Solve NLP
  status = KTR_solve(kc_handle_,
                   &output(NLP_X_OPT).front(),
                   &lambda[0],
                   0,  // not used
                   &output(NLP_COST).front(),
                   0,  // not used
                   0,  // not used
                   0,  // not used
                   0,  // not used
                   0,  // not used
                   this); // to be retrieved in the callback function
  casadi_assert(status<=0); // make sure the NLP finished solving
    
  // Copy lagrange multipliers
  output(NLP_LAMBDA_OPT).set(&lambda[0]);
  output(NLP_LAMBDA_LBX).set(&lambda[m_]);
  output(NLP_LAMBDA_UBX).set(&lambda[m_]);

  // Free memory (move to destructor!)
  KTR_free(&kc_handle_);
  kc_handle_ = 0;

}


int KnitroInternal::callback(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH, const double* const x,
                             const double* const lambda, double* const obj, double* const c, double* const objGrad,
                             double* const jac, double* const hessian, double* const hessVector, void *userParams){
  try{
    // Get a pointer to the calling object
    KnitroInternal* this_ = static_cast<KnitroInternal*>(userParams);
    
    // Direct to the correct function
    switch(evalRequestCode){
      case KTR_RC_EVALFC: this_->evalfc(x,*obj,c); break;
      case KTR_RC_EVALGA: this_->evalga(x,objGrad,jac); break;
      case KTR_RC_EVALH:  this_->evalh(x,lambda,hessian); break;
      default: casadi_assert_message(0,"KnitroInternal::callback: unknown method");
    }
    
    return 0;
  } catch (exception& ex){
    cerr << "KnitroInternal::callback caugth exception: " << ex.what() << endl;
    return -1;
  }
}

void KnitroInternal::evalfc(const double* x, double& obj, double *c){
  // Pass the argument to the function
  F_.setInput(x);

  // Evaluate the function
  F_.evaluate();

  // Get the result
  F_.getOutput(obj);
  
  // Pass the argument to the function
  G_.setInput(x);

  // Evaluate the function
  G_.evaluate();

  // Get the result
  G_.getOutput(c);
}

void KnitroInternal::evalga(const double* x, double* objGrad, double* jac){
  // Pass the argument and adjoint seed to the function
  F_.setInput(x);
  F_.setAdjSeed(1.0);

  // Evaluate the function using adjoint mode AD
  F_.evaluate(0,1);

  // Get the result
  F_.getAdjSens(objGrad);
  
  // Pass the argument to the Jacobian function
  J_.setInput(x);
  
  // Evaluate the Jacobian function
  J_.evaluate();
  
  // Get the result
  J_.getOutput(jac);
}

void KnitroInternal::evalh(const double* x, const double* lambda, double* hessian){
  // Pass input
  H_.setInput(x);
  H_.setInput(lambda,1);
  H_.setInput(1.0,2);

  // Evaluate
  H_.evaluate();

  // Get results
  H_.output().get(hessian,SPARSESYM);
}



} // namespace CasADi
