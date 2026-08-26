/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_NLPSOL_HPP
#define CASADI_NLPSOL_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_nlpsol Title
      \par

      Create an NLP solver
      Creates a solver for the following parametric nonlinear program (NLP):
      \verbatim

      min          F(x, p) + F_s(s, p)
      x, s

      subject to
      LBX - S_x s_l <=   x    <= UBX + S_x s_u
      LBG - S_g s_l <= G(x, p) <= UBG + S_g s_u
      0             <=   s    <= UBS
      p  == P

      with s = [s_l; s_u] and S = [S_g; S_x]

      nx: number of decision variables
      ng: number of constraints
      np: number of parameters
      ns: number of slacks (per side); s is 2*ns x 1

      \endverbatim

      The slack part is optional: without the 'S' option, ns is zero and the
      formulation reduces to the classical NLP.

      The slack variables themselves are declared through the 'nlp' dictionary:
      \verbatim
      s   : slack variables, 2*ns x 1, holding [s_l; s_u]
      f_s : slack objective, scalar, may depend on s and p only
      \endverbatim

      Which rows they relax is declared through the 'S' option, a Sparsity of
      shape (ng+nx)-by-ns holding [S_g; S_x]. It is structural only -- its
      entries count as 1 -- and a row of S that is structurally empty leaves
      the corresponding constraint / simple bound hard.

      The solver reports the total objective F + F_s in 'f'; 'g' and 'lam_g'
      refer to the constraints as written, not to the relaxed ones, and the
      multipliers of the slack bounds are reported in 'lam_s'.

      By default the slack layer is de-sugared into an augmented plain NLP
      before the problem reaches the plugin. A plugin that implements the
      soft-constraint formulation itself gets the structure handed over intact
      instead; the 'expand_slacks' option overrides that choice either way.
      The user-facing behaviour is identical in both cases.

      \generalsection{Nlpsol}
      \pluginssection{Nlpsol}

      \author Joel Andersson
      \date 2011-2015

      \identifier{21q} */

  /** \defgroup nlpsol Title
  * @copydoc main_nlpsol
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_nlpsol
  * \endif
  */
  ///@{
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const SXDict& nlp, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const MXDict& nlp, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const std::string& fname, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const Importer& compiler, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const NlpBuilder& nl, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const Function& nlp, const Dict& opts=Dict());
  ///@}

  /** \brief Get input scheme of NLP solvers

  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif

      \identifier{1sy} */
  CASADI_EXPORT std::vector<std::string> nlpsol_in();

  /** \brief Get NLP solver output scheme of NLP solvers

  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif

      \identifier{1sz} */
  CASADI_EXPORT std::vector<std::string> nlpsol_out();

  /** \brief Get NLP solver input scheme name by index

  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif

      \identifier{1t0} */
  CASADI_EXPORT std::string nlpsol_in(casadi_int ind);

  /** \brief Get output scheme name by index

  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif

      \identifier{1t1} */
  CASADI_EXPORT std::string nlpsol_out(casadi_int ind);

  /** \brief Number of NLP solver inputs

      \identifier{1t2} */
  CASADI_EXPORT casadi_int nlpsol_n_in();

  /** \brief Number of NLP solver outputs

      \identifier{1t3} */
  CASADI_EXPORT casadi_int nlpsol_n_out();

  ///@{
  /** \brief Default input for an NLP solver

      \identifier{1t4} */
  CASADI_EXPORT double nlpsol_default_in(casadi_int ind);
  CASADI_EXPORT std::vector<double> nlpsol_default_in();
  ///@}

  /** \brief Get all options for a plugin

      \identifier{1t5} */
  CASADI_EXPORT std::vector<std::string> nlpsol_options(const std::string& name);

  /** \brief Get type info for a particular option

      \identifier{1t6} */
  CASADI_EXPORT std::string nlpsol_option_type(const std::string& name, const std::string& op);

  /** \brief Get documentation for a particular option

      \identifier{1t7} */
  CASADI_EXPORT std::string nlpsol_option_info(const std::string& name, const std::string& op);

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_nlpsol(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_nlpsol(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_nlpsol(const std::string& name);

  /** @} */

#ifndef SWIG
/// Input arguments of an NLP function
enum NLPInput {
  /// Decision variable
  NL_X,
  /// Fixed parameter
  NL_P,
  /// Number of NLP inputs
  NL_NUM_IN
};

/// Shortname for onput arguments of an NLP function
const std::vector<std::string> NL_INPUTS = {"x", "p"};

/// Output arguments of an NLP function
enum NLPOutput {
  /// Objective function
  NL_F,
  /// Constraint function
  NL_G,
  /// Number of NLP outputs
  NL_NUM_OUT
};

/// Shortname for output arguments of an NLP function
const std::vector<std::string> NL_OUTPUTS = {"f", "g"};

/** \brief Extra IO of an NLP function in native-slack mode

    A plugin that declares it handles the slack ("soft constraint") layer
    natively is handed an oracle with one extra input and one extra output:

      nlp : (x, p, s) -> (f, g, f_s)

    'x', 'f' and 'g' are then byte-identical to the hard problem -- 'f' is the
    user's objective WITHOUT the slack penalty, which lives in 'f_s'.

    These constants deliberately live next to, rather than inside, NLPInput /
    NLPOutput: NL_NUM_IN / NL_NUM_OUT are used to size call vectors elsewhere
    (see Conic::create_rqp) and must keep meaning "plain NLP".

    \identifier{2k0} */
///@{
/// Slack variables [s_l;s_u] (2*ns x 1), input index of an NLP function
const casadi_int NL_S = 2;
/// Slack objective (scalar), output index of an NLP function
const casadi_int NL_F_S = 2;
/// Number of NLP inputs/outputs in native-slack mode
const casadi_int NL_NUM_IN_S = 3;
const casadi_int NL_NUM_OUT_S = 3;
/// Shortname for input arguments of an NLP function in native-slack mode
const std::vector<std::string> NL_INPUTS_S = {"x", "p", "s"};
/// Shortname for output arguments of an NLP function in native-slack mode
const std::vector<std::string> NL_OUTPUTS_S = {"f", "g", "f_s"};
///@}

/// Input arguments of an NLP Solver
enum NlpsolInput {
  /// Decision variables, initial guess (nx x 1)
  NLPSOL_X0,
  /// Value of fixed parameters (np x 1)
  NLPSOL_P,
  /// Decision variables lower bound (nx x 1), default -inf
  NLPSOL_LBX,
  /// Decision variables upper bound (nx x 1), default +inf
  NLPSOL_UBX,
  /// Constraints lower bound (ng x 1), default -inf
  NLPSOL_LBG,
  /// Constraints upper bound (ng x 1), default +inf
  NLPSOL_UBG,
  /// Lagrange multipliers for bounds on X, initial guess (nx x 1)
  NLPSOL_LAM_X0,
  /// Lagrange multipliers for bounds on G, initial guess (ng x 1)
  NLPSOL_LAM_G0,
  /// Slack variables, initial guess (2*ns x 1)
  NLPSOL_S0,
  /// Slack variables upper bound (2*ns x 1), default +inf
  NLPSOL_UBS,
  /// Lagrange multipliers for bounds on S, initial guess (2*ns x 1)
  NLPSOL_LAM_S0,
  NLPSOL_NUM_IN
};

/// Output arguments of an NLP Solver
enum NlpsolOutput {
  /// Decision variables at the optimal solution (nx x 1)
  NLPSOL_X,
  /// Cost function value at the optimal solution (1 x 1)
  NLPSOL_F,
  /// Constraints function at the optimal solution (ng x 1)
  NLPSOL_G,
  /// Lagrange multipliers for bounds on X at the solution (nx x 1)
  NLPSOL_LAM_X,
  /// Lagrange multipliers for bounds on G at the solution (ng x 1)
  NLPSOL_LAM_G,
  /// Lagrange multipliers for bounds on P at the solution (np x 1)
  NLPSOL_LAM_P,
  /// Slack variables at the optimal solution (2*ns x 1)
  NLPSOL_S,
  /// Lagrange multipliers for bounds on S at the solution (2*ns x 1)
  NLPSOL_LAM_S,
  NLPSOL_NUM_OUT
};
#endif // SWIG

} // namespace casadi

#endif // CASADI_NLPSOL_HPP
