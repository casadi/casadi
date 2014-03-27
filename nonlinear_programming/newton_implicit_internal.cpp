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
#include "symbolic/function/mx_function.hpp"

#include "symbolic/profiling.hpp"
#include "symbolic/casadi_options.hpp"

using namespace std;
namespace CasADi {

  NewtonImplicitInternal::NewtonImplicitInternal(const Function& f, const Function& jac, const LinearSolver& linsol) : ImplicitFunctionInternal(f,jac,linsol) {
    addOption("abstol",                      OT_REAL,1e-12,"Stopping criterion tolerance on max(|F|)");
    addOption("abstolStep",                  OT_REAL,1e-12,"Stopping criterion tolerance on step size");
    addOption("max_iter",  OT_INTEGER, 1000, "Maximum number of Newton iterations to perform before returning.");
    addOption("monitor",   OT_STRINGVECTOR, GenericType(),  "", "step|stepsize|J|F|normF", true);
  }

  NewtonImplicitInternal::~NewtonImplicitInternal(){ 
  }

  void NewtonImplicitInternal::solveNonLinear() {
    casadi_log("NewtonImplicitInternal::solveNonLinear:begin");
    
    // Set up timers for profiling
    double time_zero;
    double time_start;
    double time_stop;
    if (CasadiOptions::profiling) {
      time_zero = getRealTime();
      CasadiOptions::profilingLog  << "start " << this << ":" <<getOption("name") << std::endl; 
    }
    
    // Pass the inputs to J
    for(int i=0; i<getNumInputs(); ++i){
      if(i!=iin_) jac_.setInput(input(i),i);
    }

    // Aliases
    DMatrix &u = output(iout_);
    DMatrix &J = jac_.output(0);
    DMatrix &F = jac_.output(1+iout_);
  
    // Perform the Newton iterations
    int iter=0;
    
    bool success = true;
    
    while(true){
      // Break if maximum number of iterations already reached
      if (iter >= max_iter_) {
        log("evaluate","Max. iterations reached.");
        stats_["return_status"] = "max_iteration_reached";
        success = false;
        break;
      }

      // Start a new iteration
      iter++;

      // Print progress
      if (monitored("step") || monitored("stepsize")) {
        std::cout << "Step " << iter << "." << std::endl; 
      }
    
      if (monitored("step")) {
        std::cout << "  u = " << u << std::endl;
      }
    
      // Use u to evaluate J
      jac_.setInput(u,iin_);
      for(int i=0; i<getNumInputs(); ++i)
        if(i!=iin_) jac_.setInput(input(i),i);
      
      if (CasadiOptions::profiling) {
        time_start = getRealTime(); // Start timer
      }
    
      jac_.evaluate();
      
      // Write out profiling information
      if (CasadiOptions::profiling) {
        time_stop = getRealTime(); // Stop timer
        CasadiOptions::profilingLog  << double(time_stop-time_start)*1e6 << " ns | " << double(time_stop-time_zero)*1e3 << " ms | " << this << ":" << getOption("name") << ":0|" << jac_.get() << ":" << jac_.getOption("name") << "|evaluate jacobian" << std::endl;
      }
    
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
      linsol_.setInput(J,LINSOL_A);
      
      if (CasadiOptions::profiling) {
        time_start = getRealTime(); // Start timer
      }
      linsol_.prepare();
      // Write out profiling information
      if (CasadiOptions::profiling) {
        time_stop = getRealTime(); // Stop timer
        CasadiOptions::profilingLog  << double(time_stop-time_start)*1e6 << " ns | " << double(time_stop-time_zero)*1e3 << " ms | " << this << ":" << getOption("name") << ":1||prepare linear system" << std::endl;
      }

      if (CasadiOptions::profiling) {
        time_start = getRealTime(); // Start timer
      }
      // Solve against F
      linsol_.solve(&F.front(),1,false);
      if (CasadiOptions::profiling) {
        time_stop = getRealTime(); // Stop timer
        CasadiOptions::profilingLog  << double(time_stop-time_start)*1e6 << " ns | " << double(time_stop-time_zero)*1e3 << " ms | " << this << ":" << getOption("name") << ":2||solve linear system" << std::endl;
      }
      
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
      std::transform(u.begin(), u.end(), F.begin(), u.begin(), std::minus<double>());

      // Get auxiliary outputs
      for(int i=0; i<getNumOutputs(); ++i){
        if(i!=iout_) jac_.getOutput(output(i),1+i);
      }
    }
  
    // Store the iteration count
    if (gather_stats_) stats_["iter"] = iter; 
    
    if (success) stats_["return_status"] = "success";
  
    // Factorization up-to-date
    fact_up_to_date_ = true;
    
    casadi_log("NewtonImplicitInternal::solveNonLinear():end after " << iter << " steps");
  }

  void NewtonImplicitInternal::init(){
    
    // Call the base class initializer
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

