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

#ifndef NLP_SOLVER_INTERNAL_HPP
#define NLP_SOLVER_INTERNAL_HPP

#include "nlp_solver.hpp"
#include "fx_internal.hpp"

namespace CasADi{
    
/** \brief NLP solver storage class

  @copydoc NLPSolver_doc
  \author Joel Andersson 
  \date 2010-2013
*/
  class NLPSolverInternal : public FXInternal{

  public:
    /// Constructor
    NLPSolverInternal(const FX& nlp);

    /// Destructor
    virtual ~NLPSolverInternal() = 0;

    /// Initialize
    virtual void init();

    /// Prints out a human readable report about possible constraint violations - all constraints
    void reportConstraints(std::ostream &stream=std::cout);
  
    /// Warns the user about inital bounds, if option 'warn_initial_bounds' is true
    virtual void checkInitialBounds();
  
    /// Set options that make the NLP solver more suitable for solving QPs
    virtual void setQPOptions() { };

    /// Get or generate a function to calculate the gradient of the objective function
    virtual FX getGradF();
  
    /// Get or generate a function to calculate the Jacobian of the constraint function
    virtual FX getJacG();

    /// Get or generate a function to calculate the gradient of the Lagrangian function
    virtual FX getGradLag();

    /// Get or generate a function to calculate the Hessian of the Lagrangian function
    virtual FX getHessLag();

    /// Get or generate the sparsity pattern of the Hessian of the Lagrangian
    virtual CRSSparsity getHessLagSparsity();
    
    // Access the objective gradient function
    FX& gradF();

    /// Access the Jacobian of the constraint function
    FX& jacG();

    /// Access the Hessian of the Lagrangian function
    FX& hessLag();

    /// Access the gradient of the Lagrangian function
    FX& gradLag();

    /// Get the sparsity pattern of the Hessian of the Lagrangian
    CRSSparsity& hessLagSparsity();

    /// Number of variables
    int nx_;
  
    /// Number of constraints
    int ng_;
  
    /// Number of parameters
    int np_;
  
    /// callback function, executed at each iteration
    FX callback_;
  
    /// Execute the callback function only after this amount of iterations
    int callback_step_;
  
    /// The NLP
    FX nlp_;

    // Gradient of the objective
    FX gradF_;
    
    // Jacobian of the constraints
    FX jacG_;
    
    // Hessian of the Lagrangian
    FX hessLag_;

    // Gradient of the Lagrangian
    FX gradLag_;

    // Sparsity pattern of the Hessian of the Lagrangian
    CRSSparsity hessLagSparsity_;
  };

} // namespace CasADi

#endif //NLP_SOLVER_INTERNAL_HPP
