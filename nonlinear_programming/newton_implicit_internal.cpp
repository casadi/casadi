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
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/fx/mx_function.hpp"

using namespace std;
namespace CasADi {

  NewtonImplicitInternal::NewtonImplicitInternal(const FX& f, const FX& jac, const LinearSolver& linsol) : ImplicitFunctionInternal(f,jac,linsol) {
    addOption("abstol",                      OT_REAL,1e-12,"Stopping criterion tolerance on max(|F|)");
    addOption("abstolStep",                  OT_REAL,1e-12,"Stopping criterion tolerance on step size");
    addOption("max_iter",  OT_INTEGER, 1000, "Maximum number of Newton iterations to perform before returning.");
    addOption("monitor",   OT_STRINGVECTOR, GenericType(),  "", "step|stepsize|J|F|normF", true);
  }

  NewtonImplicitInternal::~NewtonImplicitInternal(){ 
  }

  double Xk_update (double Xk, double step) { return Xk-step; }

  void NewtonImplicitInternal::evaluate(int nfdir, int nadir) {
    casadi_log("NewtonImplicitInternal::evaluate(" << nfdir << ", " << nadir<< "):begin");
    // Pass the inputs to J
    for (int i=1;i<jac_.getNumInputs();++i) {
      std::copy(input(i-1).data().begin(),input(i-1).data().end(),jac_.input(i).data().begin());
    }

    // Aliases
    DMatrix &Xk = output(0);
    DMatrix &J = jac_.output(0);
    DMatrix &F = jac_.output(1);
  
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
      std::copy(Xk.data().begin(),Xk.data().end(),jac_.input().data().begin());
      jac_.evaluate();
    
      if (monitored("F")) std::cout << "  F = " << F << std::endl;
      if (monitored("normF")) std::cout << "  F (min, max, 1-norm, 2-norm) = " << (*std::min_element(F.data().begin(),F.data().end())) << ", " << (*std::max_element(F.data().begin(),F.data().end())) << ", " << sumAll(fabs(F)) << ", " << sqrt(sumAll(F*F)) << std::endl;
      if (monitored("J")) std::cout << "  J = " << J << std::endl;

      if ( numeric_limits<double>::infinity() != abstol_ ) {
	double maxF = std::max((*std::max_element(F.data().begin(),F.data().end())),-(*std::min_element(F.data().begin(),F.data().end())));
	if (maxF <= abstol_) {
	  casadi_log("Converged to acceptable tolerance - abstol: " << abstol_);
	  break;
	}
      } 
    
      // Prepare the linear solver with J
      linsol_.setInput(J,0);
      linsol_.prepare();
    
      // Solve against F
      linsol_.solve(&F.front(),1,false);

      if (monitored("step")) {
	std::cout << "  step = " << F << std::endl;
      }
    
      if ( numeric_limits<double>::infinity() != abstolStep_ ) {
	double maxF = std::max((*std::max_element(F.data().begin(),F.data().end())),-(*std::min_element(F.data().begin(),F.data().end())));
	if (monitored("stepsize")) {
	  std::cout << "  stepsize = " << maxF << std::endl;
	}
	if (maxF <= abstolStep_) {
	  casadi_log("Converged to acceptable tolerance - abstolStep: " << abstolStep_);
	  break;
	}
      } 
    
      // Update Xk+1 = Xk - J^(-1) F
      std::transform(Xk.begin(), Xk.end(), F.begin(), Xk.begin(), Xk_update);
  
    }
  
    int Niter = i+1;
    if (gather_stats_) stats_["iter"] = Niter; 
  
    // Pass the remainder of outputs
    for (int i=2;i<jac_.getNumOutputs();++i) {
      std::copy(jac_.output(i).data().begin(),jac_.output(i).data().end(),output(i-1).data().begin());
    }
  
    // Delegate calculation of sensitivities to base class
    if(nfdir!=0 || nadir!=0)
      evaluate_sens(nfdir,nadir,true);
  
    casadi_log("NewtonImplicitInternal::evaluate(" << nfdir << ", " << nadir<< "):end after " << Niter << " steps");
  }

  void NewtonImplicitInternal::init(){

    ImplicitFunctionInternal::init();

    casadi_assert_message(f_.getNumInputs()>0,"NewtonImplicitInternal: the supplied f must have at least one input.");
  
    casadi_assert_message(!linsol_.isNull(),"NewtonImplicitInternal::init: linear_solver must be supplied");
  
    if (hasSetOption("max_iter"))
      max_iter_ = getOption("max_iter");
    
    if (hasSetOption("abstol"))
      abstol_ = getOption("abstol");

    if (hasSetOption("abstolStep"))
      abstolStep_ = getOption("abstolStep");
    
  }

} // namespace CasADi

