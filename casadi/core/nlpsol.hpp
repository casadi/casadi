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


#ifndef CASADI_NLPSOL_HPP
#define CASADI_NLPSOL_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_nlpsol
      Create an NLP solver
      Creates a solver for the following parametric nonlinear program (NLP):
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

      \generalsection{Nlpsol}
      \pluginssection{Nlpsol}

      \author Joel Andersson
      \date 2011-2015
  */

  /** \defgroup nlpsol
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
#ifndef SWIG
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const Function& nlp, const Dict& opts=Dict());
#endif // SWIG
  ///@}

  /** \brief Get input scheme of NLP solvers
  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif
  */
  CASADI_EXPORT std::vector<std::string> nlpsol_in();

  /** \brief Get NLP solver output scheme of NLP solvers
  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif
  */
  CASADI_EXPORT std::vector<std::string> nlpsol_out();

  /** \brief Get NLP solver input scheme name by index
  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif
  */
  CASADI_EXPORT std::string nlpsol_in(int ind);

  /** \brief Get output scheme name by index
  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif
  */
  CASADI_EXPORT std::string nlpsol_out(int ind);

  /** \brief Number of NLP solver inputs */
  CASADI_EXPORT int nlpsol_n_in();

  /** \brief Number of NLP solver outputs */
  CASADI_EXPORT int nlpsol_n_out();

  ///@{
  /** \brief Default input for an NLP solver */
  CASADI_EXPORT double nlpsol_default_in(int ind);
  CASADI_EXPORT std::vector<double> nlpsol_default_in();
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_nlpsol(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_nlpsol(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_nlpsol(const std::string& name);

  /** @} */

} // namespace casadi

#endif // CASADI_NLPSOL_HPP
