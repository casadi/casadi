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

#include "newton_implicit_internal.hpp"

#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/fx/mx_function.hpp"

using namespace std;
namespace CasADi {

NewtonImplicitInternal* NewtonImplicitInternal::clone() const{
  // Return a deep copy
  NewtonImplicitInternal* node = new NewtonImplicitInternal(f_,nrhs_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
NewtonImplicitInternal::NewtonImplicitInternal(const FX& f, int nrhs) : ImplicitFunctionInternal(f,nrhs) {
  addOption("abstol",                      OT_REAL,1e-12,"Stopping criterion tolerance");
  addOption("max_iter",  OT_INTEGER, 1000, "Maximum number of Newton iterations to perform before returning.");
  addOption("monitor",   OT_STRINGVECTOR, GenericType(),  "", "step|stepsize|J|F", true);
}

NewtonImplicitInternal::~NewtonImplicitInternal(){ 
}

double Xk_update (double Xk, double step) { return Xk-step; }

void NewtonImplicitInternal::evaluate(int nfdir, int nadir) {

  // Pass the inputs to J
  for (int i=1;i<J_.getNumInputs();++i) {
    std::copy(input(i-1).data().begin(),input(i-1).data().end(),J_.input(i).data().begin());
  }

  // Aliases
  DMatrix &Xk = output(0);
  DMatrix &J = J_.output(0);
  DMatrix &F = J_.output(1);
  
  int i=0;
  for (;true;++i) {
    if (i>=max_iter_) {
      log("evaluate","Max. iterations reached.");
      break;
    }
    if (monitored("step") || monitored("stepsize")) {
      std::cout << "Step " << i << "." << std::endl; 
    }
    
    if (monitored("step")) {
      std::cout << "  Xk = " << Xk << std::endl;
    }
    
    // Use Xk to evaluate J
    std::copy(Xk.data().begin(),Xk.data().end(),J_.input().data().begin());
    J_.evaluate();
    
    if (monitored("F")) std::cout << "  F = " << F << std::endl;
    if (monitored("J")) std::cout << "  J = " << J << std::endl;
    
    // Prepare the linear solver with J
    linsol_.setInput(J,0);
    linsol_.prepare();
    
    // Solve against F
    linsol_.solve(&F.front(),1,false);

    if (monitored("step")) {
      std::cout << "  step = " << F << std::endl;
    }
    
    if ( numeric_limits<double>::infinity() != abstol_ ) {
      double maxF = 0;
      for (int k=0;k<F.size();++k) {
        if (F.data()[k]>maxF)  maxF = F.data()[k];
        if (-F.data()[k]>maxF) maxF = -F.data()[k];
      }
      if (monitored("stepsize")) {
        std::cout << "  stepsize = " << maxF << std::endl;
      }
      if (maxF <= abstol_) {
        log("evaluate","Converged to acceptable tolerance");
        break;
      }
    } 
    
    // Update Xk+1 = Xk - J^(-1) F
    std::transform(Xk.begin(), Xk.end(), F.begin(), Xk.begin(), Xk_update);
  
  }
  
  if (gather_stats_) stats_["iter"] = i+1; 
  
  // Pass the remainder of outputs
  for (int i=2;i<J_.getNumOutputs();++i) {
    std::copy(J_.output(i).data().begin(),J_.output(i).data().end(),output(i-1).data().begin());
  }
  
  // End of function if no sensitivities
  if(nfdir==0 && nadir==0)
    return;
  
  evaluate_sens(nfdir,nadir,true);

}

void NewtonImplicitInternal::init(){

  ImplicitFunctionInternal::init();

  casadi_assert_message(f_.getNumInputs()>0,"NewtonImplicitInternal: the supplied f must have at least one input.");
  
  casadi_assert_message(!linsol_.isNull(),"NewtonImplicitInternal::init: linear_solver must be supplied");
  
  if (hasSetOption("max_iter"))
    max_iter_ = getOption("max_iter");
    
  if (hasSetOption("abstol"))
    abstol_ = getOption("abstol");
}

} // namespace CasADi

