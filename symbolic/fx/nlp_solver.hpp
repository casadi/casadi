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

#ifndef NLP_SOLVER_HPP
#define NLP_SOLVER_HPP

#include "fx.hpp"


/** \defgroup NLPSolver_doc

  Solves the following parametric nonlinear program (NLP):
  \verbatim
  min          F(x,p)
   x
  
  subject to
              LBX <=   x    <= UBX
              LBG <= G(x,p) <= UBG
                         p  == P
              
      nx: number of decision variables
      ng: number of constraints
      np: number of parameters
  \endverbatim

*/

namespace CasADi{

  /// Input arguments of an NLP function [nlpIn]
  enum NLPInput{
    /// Decision variable [x]
    NLP_X,
    /// Fixed parameter [p]
    NLP_P, 
    /// Number of NLP inputs
    NLP_NUM_IN
  };

  /// Output arguments of an NLP function [nlpOut]
  enum NLPOutput{ 
    /// Objective function [f]
    NLP_F,
    /// Constraint function [g]
    NLP_G,
    /// Number of NLP outputs
    NLP_NUM_OUT
  };

  /// Input arguments of an NLP objective gradient function [gradFIn]
  enum GradFInput{
    /// Decision variable [x]
    GRADF_X,
    /// Fixed parameter [p]
    GRADF_P, 
    /// Number of inputs
    GRADF_NUM_IN
  };

  /// Output arguments of an NLP objective gradient function [gradFOut]
  enum GradFOutput{ 
    /// Jacobian of the constraints [grad]
    GRADF_GRAD,
    /// Objective function [f]
    GRADF_F,
    /// Constraint function [g]
    GRADF_G,
    /// Number of outputs
    GRADF_NUM_OUT
  };

  /// Input arguments of an NLP Jacobian function [jacGIn]
  enum JacGInput{
    /// Decision variable [x]
    JACG_X,
    /// Fixed parameter [p]
    JACG_P, 
    /// Number of inputs
    JACG_NUM_IN
  };

  /// Output arguments of an NLP Jacobian function [jacGOut]
  enum JacGOutput{ 
    /// Jacobian of the constraints [jac]
    JACG_JAC,
    /// Objective function [f]
    JACG_F,
    /// Constraint function [g]
    JACG_G,
    /// Number of outputs
    JACG_NUM_OUT
  };

  /// Input arguments of an NLP Hessian function [hessLagIn]
  enum HessLagInput{
    /// Decision variable [x]
    HESSLAG_X,
    /// Fixed parameter [p]
    HESSLAG_P, 
    /// Multiplier for f [lam_f]
    HESSLAG_LAM_F,
    /// Multiplier for g [lam_g]
    HESSLAG_LAM_G,
    /// Number of inputs
    HESSLAG_NUM_IN
  };

  /// Output arguments of an NLP Hessian function [hessLagOut]
  enum HessLagOutput{ 
    /// Hessian of the Lagrangian [hess]
    HESSLAG_HESS,
    /// Objective function [f]
    HESSLAG_F,
    /// Constraint function [g]
    HESSLAG_G,
    /// Gradient of the Lagrangian with respect to x [grad_x]
    HESSLAG_GRAD_X,
    /// Gradient of the Lagrangian with respect to p [grad_p]
    HESSLAG_GRAD_P,
    /// Number of outputs
    HESSLAG_NUM_OUT
  };
  
  /// Input arguments of an NLP Solver [nlpSolverIn]
  enum NLPSolverInput{
    /// Decision variables, initial guess (nx x 1)  [x0]
    NLP_SOLVER_X0,
    /// Value of fixed parameters (np x 1) [p]
    NLP_SOLVER_P,
    /// Decision variables lower bound (nx x 1), default -inf [lbx]
    NLP_SOLVER_LBX,
    /// Decision variables upper bound (nx x 1), default +inf [ubx]
    NLP_SOLVER_UBX,
    /// Constraints lower bound (ng x 1), default -inf [lbg]
    NLP_SOLVER_LBG,
    /// Constraints upper bound (ng x 1), default +inf [ubg]
    NLP_SOLVER_UBG,
    /// Lagrange multipliers for bounds on X, initial guess (nx x 1) [lam_x0]
    NLP_SOLVER_LAM_X0,
    /// Lagrange multipliers for bounds on G, initial guess (ng x 1) [lam_g0]
    NLP_SOLVER_LAM_G0,
    NLP_SOLVER_NUM_IN
  };

  /// Output arguments of an NLP Solver [nlpSolverOut]
  enum NLPSolverOutput{
    /// Decision variables at the optimal solution (nx x 1) [x]
    NLP_SOLVER_X,
    /// Cost function value at the optimal solution (1 x 1) [f]
    NLP_SOLVER_F,
    /// Constraints function at the optimal solution (ng x 1) [g]
    NLP_SOLVER_G,
    /// Lagrange multipliers for bounds on X at the solution (nx x 1) [lam_x]
    NLP_SOLVER_LAM_X, 
    /// Lagrange multipliers for bounds on G at the solution (ng x 1) [lam_g]
    NLP_SOLVER_LAM_G,
    /// Lagrange multipliers for bounds on P at the solution (np x 1) [lam_p]
    NLP_SOLVER_LAM_P, 
    NLP_SOLVER_NUM_OUT
  };

  class NLPSolverInternal;

  /** \brief NLPSolver

      @copydoc NLPSolver_doc

      \author Joel Andersson 
      \date 2010
  */
  class NLPSolver : public FX{
  public:

    /// Default constructor
    NLPSolver();

    /// Access functions of the node
    NLPSolverInternal* operator->();
    const NLPSolverInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  
    /// Prints out a human readable report about possible constraint violations, after solving 
    void reportConstraints(std::ostream &stream=std::cout);

    std::string getReportConstraints() { std::stringstream s; reportConstraints(s); return s.str(); }
  
    /// Set options that make the NLP solver more suitable for solving QPs
    void setQPOptions();
  
    /// Access the NLP
    FX nlp();

    // Access the objective gradient function
    FX gradF();

    /// Access the Jacobian of the constraint function
    FX jacG();

    /// Access the Hessian of the Lagrangian function
    FX hessLag();

    /// Join F and G in old signature style to a common NLP function
    static FX joinFG(FX F, FX G);
  };

} // namespace CasADi

#endif // NLP_SOLVER_HPP

