/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_NLP_SOLVER_HPP
#define CASADI_NLP_SOLVER_HPP

#include "function.hpp"


/** \defgroup NlpSolver_doc

  Solves the following parametric nonlinear program (NLP):
  \verbatim
  min          F(x, p)
   x

  subject to
              LBX <=   x    <= UBX
              LBG <= G(x, p) <= UBG
                         p  == P

      nx: number of decision variables
      ng: number of constraints
      np: number of parameters
  \endverbatim

*/

namespace casadi {

  /// Input arguments of an NLP function [nlpIn]
  enum NLPInput {
    /// Decision variable [x]
    NL_X,
    /// Fixed parameter [p]
    NL_P,
    /// Number of NLP inputs
    NL_NUM_IN
  };

  /// Output arguments of an NLP function [nlpOut]
  enum NLPOutput {
    /// Objective function [f]
    NL_F,
    /// Constraint function [g]
    NL_G,
    /// Number of NLP outputs
    NL_NUM_OUT
  };

  /// Input arguments of an NLP objective gradient function [gradFIn]
  enum GradFInput {
    /// Decision variable [x]
    GRADF_X,
    /// Fixed parameter [p]
    GRADF_P,
    /// Number of inputs
    GRADF_NUM_IN
  };

  /// Output arguments of an NLP objective gradient function [gradFOut]
  enum GradFOutput {
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
  enum JacGInput {
    /// Decision variable [x]
    JACG_X,
    /// Fixed parameter [p]
    JACG_P,
    /// Number of inputs
    JACG_NUM_IN
  };

  /// Output arguments of an NLP Jacobian function [jacGOut]
  enum JacGOutput {
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
  enum HessLagInput {
    /// Decision variable [x]
    HESSLAG_X,
    /// Fixed parameter [p]
    HESSLAG_P,
    /// Multiplier for f. Just a scalar factor for the objective that the
    /// NLP solver might use to scale the objective. [lam_f]
    HESSLAG_LAM_F,
    /// Multiplier for g [lam_g]
    HESSLAG_LAM_G,
    /// Number of inputs
    HESSLAG_NUM_IN
  };

  /// Output arguments of an NLP Hessian function [hessLagOut]
  enum HessLagOutput {
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
  enum NlpSolverInput {
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
  enum NlpSolverOutput {
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

  class NlpSolverInternal;

  /** \brief NlpSolver

      @copydoc NlpSolver_doc

      \generalsection{NlpSolver}
      \pluginssection{NlpSolver}

      \author Joel Andersson
      \date 2010
  */
  class CASADI_CORE_EXPORT NlpSolver : public Function {
  public:

    /// Default constructor
    NlpSolver();

    /// NLP solver factory
    NlpSolver(
      const std::string& name,
      /**< \pluginargument{NlpSolver}
      */
      const Function& nlp
      /**< \parblock
       *  nlp function: \f$ [\mathbb {R}^{n_x} \times \mathbb{R}^{n_p}]
       * \mapsto [\mathbb {R} \times \mathbb{R}^{n_g}]\f$
       *
       *  @copydoc scheme_NLPInput
       *  @copydoc scheme_NLPOutput
       *
       *  \endparblock
       */
      ); // NOLINT(whitespace/parens)
    /// Access functions of the node
    NlpSolverInternal* operator->();
    const NlpSolverInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Prints out a human readable report about possible constraint violations, after solving
    void reportConstraints(std::ostream &stream=std::cout);

    std::string getReportConstraints()
    { std::stringstream s; reportConstraints(s); return s.str(); }

    /// Set options that make the NLP solver more suitable for solving QPs
    void setQPOptions();

    /** \brief Access the NLP
    *  \copydoc scheme_NlpSolverInput
    *  \copydoc scheme_NlpSolverOutput
    */
    Function nlp();

    /** Access the objective gradient function
    *  \copydoc scheme_GradFInput
    *  \copydoc scheme_GradFIOutput
    */
    Function gradF();

    /** \brief Access the Jacobian of the constraint function
    *  \copydoc scheme_HessLagInput
    *  \copydoc scheme_HessLagOutput
    */
    Function jacG();

    /** \brief Access the Hessian of the Lagrangian function
    *  \copydoc scheme_JacGInput
    *  \copydoc scheme_JacGOutput
    */
    Function hessLag();

    /// Join F and G in old signature style to a common NLP function
    static Function joinFG(Function F, Function G);

    /** \brief Get the reduced Hessian.
     * Requires a patched sIPOPT installation, see CasADi documentation. */
    DMatrix getReducedHessian();

    /// Read options from parameter xml
    void setOptionsFromFile(const std::string & file);
  };

} // namespace casadi

#endif // CASADI_NLP_SOLVER_HPP

