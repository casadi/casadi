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
  class NlpSolverInternal;

  /** \brief NlpSolver

      @copydoc NlpSolver_doc

      \generalsection{NlpSolver}
      \pluginssection{NlpSolver}

      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT NlpSolver : public Function {
  public:

    /// Default constructor
    NlpSolver();

    /// NLP solver factory (new syntax, includes initialization)
    NlpSolver(const std::string& name, const std::string& solver,
              /**< \pluginargument{NlpSolver}
               */
              const Function& nlp,
              /**< \parblock
               *  nlp function: \f$ [\mathbb {R}^{n_x} \times \mathbb{R}^{n_p}]
               * \mapsto [\mathbb {R} \times \mathbb{R}^{n_g}]\f$
               *
               *  @copydoc scheme_NLPInput
               *  @copydoc scheme_NLPOutput
               *
               *  \endparblock
               */
              const Dict& opts=Dict()
              ); // NOLINT(whitespace/parens)

    /// Access functions of the node
    NlpSolverInternal* operator->();
    const NlpSolverInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectNode* ptr);

#ifndef SWIG
    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);
#endif // SWIG

    /// Prints out a human readable report about possible constraint violations, after solving
    void reportConstraints(std::ostream &stream=casadi::userOut());

    std::string getReportConstraints()
    { std::stringstream s; reportConstraints(s); return s.str(); }

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

    /** \brief Access the Hessian of the Lagrangian function
    *  \copydoc scheme_JacGInput
    *  \copydoc scheme_JacGOutput
    */
    Function jacG();

    /** \brief Access the Jacobian of the constraint function
    *  \copydoc scheme_HessLagInput
    *  \copydoc scheme_HessLagOutput
    */
    Function hessLag();

    /** \brief Get the reduced Hessian.
     * Requires a patched sIPOPT installation, see CasADi documentation. */
    DMatrix getReducedHessian();

    /// Read options from parameter xml
    void setOptionsFromFile(const std::string & file);
  };

} // namespace casadi

#endif // CASADI_NLP_SOLVER_HPP

