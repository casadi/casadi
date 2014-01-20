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

#include "rk_base_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

  RKBaseInternal::RKBaseInternal(const FX& f, const FX& g) : IntegratorInternal(f,g){
    addOption("number_of_finite_elements",     OT_INTEGER,  20, "Number of finite elements");
    addOption("implicit_solver",               OT_IMPLICITFUNCTION,  GenericType(), "An implicit function solver");
    addOption("implicit_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLP Solver");
  }

  void RKBaseInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){    
    IntegratorInternal::deepCopyMembers(already_copied);
    implicit_solver_ = deepcopy(implicit_solver_,already_copied);
  }

  RKBaseInternal::~RKBaseInternal(){
  }

  void RKBaseInternal::init(){
    // Call the base class init
    IntegratorInternal::init();
  
    // Number of finite elements and time steps
    nk_ = getOption("number_of_finite_elements");
    casadi_assert(nk_>0);
    h_ = (tf_ - t0_)/nk_;
    
    // Setup discrete time dynamics
    setupFG();

    // Get discrete time dimensions
    Z_ = F_.input(DAE_Z);
    nZ_ = Z_.size();
    RZ_ = G_.isNull() ? DMatrix() : G_.input(RDAE_RZ);
    nRZ_ =  RZ_.size();

    // Allocate tape if backward states are present
    if(nrx_>0){
      x_tape_.resize(nk_+1,vector<double>(nx_));
      Z_tape_.resize(nk_,vector<double>(nZ_));
    }

    // Allocate a root-finding solver for the forward problem
    if(nZ_>0){

      // Get the NLP creator function
      implicitFunctionCreator implicit_function_creator = getOption("implicit_solver");
  
      // Allocate an NLP solver
      implicit_solver_ = implicit_function_creator(F_,FX(),LinearSolver());
      implicit_solver_.setOption("name",string(getOption("name")) + "_implicit_solver");
      implicit_solver_.setOption("implicit_input",DAE_Z);
      implicit_solver_.setOption("implicit_output",DAE_ALG);
    
      // Pass options
      if(hasSetOption("implicit_solver_options")){
        const Dictionary& implicit_solver_options = getOption("implicit_solver_options");
        implicit_solver_.setOption(implicit_solver_options);
      }
  
      // Initialize the solver
      implicit_solver_.init();
    }

    // Allocate a root-finding solver for the backward problem
    if(nRZ_>0){

      // Get the NLP creator function
      implicitFunctionCreator backward_implicit_function_creator = getOption("implicit_solver");
  
      // Allocate an NLP solver
      backward_implicit_solver_ = backward_implicit_function_creator(G_,FX(),LinearSolver());
      backward_implicit_solver_.setOption("name",string(getOption("name")) + "_backward_implicit_solver");
      backward_implicit_solver_.setOption("implicit_input",RDAE_RZ);
      backward_implicit_solver_.setOption("implicit_output",RDAE_ALG);
    
      // Pass options
      if(hasSetOption("implicit_solver_options")){
        const Dictionary& backward_implicit_solver_options = getOption("implicit_solver_options");
        backward_implicit_solver_.setOption(backward_implicit_solver_options);
      }
  
      // Initialize the solver
      backward_implicit_solver_.init();
    }
  }

  void RKBaseInternal::integrate(double t_out){
    // Get discrete time sought
    int k_out = std::ceil((t_out-t0_)/h_);
    k_out = std::min(k_out,nk_); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out>=0);

    // Explicit discrete time dynamics
    FX& F = nZ_>0 ? implicit_solver_ : F_;

    // Take time steps until end time has been reached
    while(k_<k_out){
      // Take step
      F.input(DAE_T).set(t_);
      F.input(DAE_X).set(output(INTEGRATOR_XF));
      F.input(DAE_Z).set(Z_);
      F.input(DAE_P).set(input(INTEGRATOR_P));
      F.evaluate();
      F.output(DAE_ODE).get(output(INTEGRATOR_XF));
      F.output(DAE_ALG).get(Z_);
      transform(F.output(DAE_QUAD).begin(),F.output(DAE_QUAD).end(),output(INTEGRATOR_QF).begin(),output(INTEGRATOR_QF).begin(),std::plus<double>());

      // Tape
      if(nrx_>0){
        output(INTEGRATOR_XF).get(x_tape_.at(k_+1));
        Z_.get(Z_tape_.at(k_));
      }

      // Advance time
      k_++;
      t_ = t0_ + k_*h_;
    }
  }

  void RKBaseInternal::integrateB(double t_out){
    // Get discrete time sought
    int k_out = std::floor((t_out-t0_)/h_);
    k_out = std::max(k_out,0); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out<=nk_);

    // Explicit discrete time dynamics
    FX& G = nRZ_>0 ? backward_implicit_solver_ : G_;

    // Take time steps until end time has been reached
    while(k_>k_out){
      // Advance time
      k_--;
      t_ = t0_ + k_*h_;

      // Take step
      G.input(RDAE_T).set(t_);
      G.input(RDAE_X).set(x_tape_.at(k_));
      G.input(RDAE_Z).set(Z_tape_.at(k_));
      G.input(RDAE_P).set(input(INTEGRATOR_P));
      G.input(RDAE_RX).set(output(INTEGRATOR_RX0));
      G.input(RDAE_RZ).set(RZ_);
      G.input(RDAE_RP).set(input(INTEGRATOR_RP));
      G.evaluate();
      G.output(RDAE_ODE).get(output(INTEGRATOR_RXF));
      G.output(RDAE_ALG).get(RZ_);
      transform(G.output(RDAE_QUAD).begin(),G.output(RDAE_QUAD).end(),output(INTEGRATOR_RQF).begin(),output(INTEGRATOR_RQF).begin(),std::plus<double>());
    }
  }

  void RKBaseInternal::reset(){
    // Reset the base classes
    IntegratorInternal::reset();

    // Bring discrete time to the beginning
    k_ = 0;

    // Get consistent initial conditions
    calculateInitialConditions();

    // Add the first element in the tape
    if(nrx_>0){
      output(INTEGRATOR_XF).get(x_tape_.at(0));
    }
  }

  void RKBaseInternal::resetB(){
    // Reset the base classes
    IntegratorInternal::resetB();

    // Bring discrete time to the end
    k_ = nk_;

    // Get consistent initial conditions
    calculateBackwardInitialConditions();
  }

  void RKBaseInternal::calculateInitialConditions(){
    casadi_assert(Z_.size()==z_.size());
    Z_.set(z_);
  }

  void RKBaseInternal::calculateBackwardInitialConditions(){
    casadi_assert(RZ_.size()==rz_.size());
    RZ_.set(rz_);
  }

} // namespace CasADi
