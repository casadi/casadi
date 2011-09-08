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

#include "nlp_solver_internal.hpp"

INPUTSCHEME(NLPInput)
OUTPUTSCHEME(NLPOutput)

using namespace std;
namespace CasADi{

NLPSolverInternal::NLPSolverInternal(const FX& F, const FX& G, const FX& H, const FX& J) : F_(F), G_(G), H_(H), J_(J){
  // set default options
  setOption("name",            "unnamed NLP solver"); // name of the function
  addOption("exact_hessian",    OT_BOOLEAN,     GenericType(), "Use an exact Hessian [default: only if one has been provided]");
  addOption("expand_f",         OT_BOOLEAN,     false,         "Expand the objective function in terms of scalar operations, i.e. MX->SX");
  addOption("expand_g",         OT_BOOLEAN,     false,         "Expand the constraint function in terms of scalar operations, i.e. MX->SX");

  n_ = 0;
  m_ = 0;
  pn_ = 0;
  pm_ = 0;
}

NLPSolverInternal::~NLPSolverInternal(){
}

void NLPSolverInternal::init(){
  // Initialize the functions
  casadi_assert_message(!F_.isNull(),"No objective function");
  if(!F_.isInit()) F_.init();
  if(!G_.isNull() && !G_.isInit()) G_.init();
  if(!H_.isNull() && !H_.isInit()) H_.init();

  // Create a Jacobian if it does not already exists
  if(!G_.isNull() && J_.isNull()){
    J_ = G_.jacobian();
  }
  if(!J_.isNull() && !J_.isInit()) J_.init();

  // Get dimensions
  n_ = F_.input(0).numel();
  m_ = G_.isNull() ? 0 : G_.output(0).numel();
    
  // Basic sanity checks
  casadi_assert_message(F_.getNumInputs()==1 or F_.getNumInputs()==2, "Wrong number of input arguments to F. Must ");
  casadi_assert_message(F_.getNumOutputs()>=1, "Wrong number of output arguments to F");
  casadi_assert_message(F_.output().scalar() && F_.output().dense(), "Output argument of F not dense scalar.");

  if(!G_.isNull()) {
    casadi_assert_message(G_.getNumInputs()>=1, "Wrong number of input arguments to G");
    casadi_assert_message(G_.getNumOutputs()>=1, "Wrong number of output arguments to G");
    casadi_assert_message(G_.input().numel()==n_, "Inconsistent dimensions");
  }
  
  if(!H_.isNull()) {
    casadi_assert_message(H_.getNumInputs()>=3, "Wrong number of input arguments to H");
    casadi_assert_message(H_.getNumOutputs()>=1, "Wrong number of output arguments to H");
    casadi_assert_message(H_.input(0).numel()==n_,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size1()==n_,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size2()==n_,"Inconsistent dimensions");
  }

  if(!J_.isNull()){
    casadi_assert_message(J_.getNumInputs()>=1, "Wrong number of input arguments to J");
    casadi_assert_message(J_.getNumOutputs()>=1, "Wrong number of output arguments to J");
    casadi_assert_message(J_.input().numel()==n_,"Inconsistent dimensions");
    casadi_assert_message(J_.output().size2()==n_,"Inconsistent dimensions");
  }
  
  pn_=0;
  pm_=0;
  
  int pn=0;
  int pm=0;
    
  // Check if any of the functions have a second argument (i.e. to pass parameters)
  FX* ff[4] = {&F_, &G_, &H_, &J_};
  for (int k=0; k<4; ++k){
    if (ff[k]->isNull()) continue;
    if (ff[k]->getNumInputs()!=2) continue;
    pn = ff[k]->input(1).size1();
    pm = ff[k]->input(1).size2();
    
    if (pn==0 or pm==0)
     continue;
    
    if ((pn!=pn_ || pm!=pm_) && pn_!=0 && pm_!=0) {
      stringstream s;
      s << "One of your supplied functions had a second input argument, which was interpreted as a parameter of shape (" << pn_ << "x" << pm_ << ")." << std::endl;
      s << "However, another function had a second input argument of shape (" << pn << "x" << pm << ")." << std::endl;
      s << "This is inconsistent." << std::endl;
      throw CasadiException(s.str());
    }
    pn_ = pn;
    pm_ = pm;
  }
  
  // Allocate space for inputs
  input_.resize(NLP_NUM_IN);
  input(NLP_X_INIT)      = DMatrix(n_,1,0);
  input(NLP_LBX)         = DMatrix(n_,1,0);
  input(NLP_UBX)         = DMatrix(n_,1,0);
  input(NLP_LBG)         = DMatrix(m_,1,0);
  input(NLP_UBG)         = DMatrix(m_,1,0);
  input(NLP_LAMBDA_INIT) = DMatrix(m_,1,0);
  input(NLP_P)           = DMatrix(pn_,pm_,0);
  
  // Allocate space for outputs
  output_.resize(NLP_NUM_OUT);
  output(NLP_X_OPT)      = DMatrix(n_,1,0);
  output(NLP_COST)       = DMatrix(1,1,0);
  output(NLP_LAMBDA_OPT) = DMatrix(m_,1,0);
  output(NLP_LAMBDA_LBX) = DMatrix(n_,1,0);
  output(NLP_LAMBDA_UBX) = DMatrix(n_,1,0);

  // Call the initialization method of the base class
  FXInternal::init();

  
  // 
//   addOption("expand_f",         OT_BOOLEAN,     false,         "Expand the objective function in terms of scalar operations, i.e. MX->SX");
//   addOption("expand_g",         OT_BOOLEAN,     false,         "Expand the constraint function in terms of scalar operations, i.e. MX->SX");
//   addOption("exact_hessian",    OT_BOOLEAN,     GenericType(), "Use an exact Hessian [default: only if one has been provided]");
//   addOption("bfgs_updates",     OT_BOOLEAN,     GenericType(), "Update the Hessian approximation with BFGS updates [default: if not exact_hessian]");

  // Read options
/*  bfgs_updates*/
  
  
}


} // namespace CasADi
