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


#ifndef CASADI_LINSOL_HPP
#define CASADI_LINSOL_HPP

#include "function.hpp"

namespace casadi {


  /** \defgroup main_linsol
   * Create a solver for linear systems of equations
   * Solves the linear system A*X = B or A^T*X = B for X
   * with A square and non-singular
   *
   *  If A is structurally singular, an error will be thrown during init.
   *  If A is numerically singular, the prepare step will fail.
   *
   *  The usual procedure to use Linsol is: \n
   *  -# init()
   *  -# set the first input (A)
   *  -# prepare()
   *  -# set the second input (b)
   *  -# solve()
   *  -# Repeat steps 4 and 5 to work with other b vectors.
   *
   * The standard evaluation combines the prepare() and solve() step and may
   * therefore more expensive if A is invariant.
   *
   *  \generalsection{Linsol}
   *  \pluginssection{Linsol}
   * \author Joel Andersson
   * \date 2011-2015
   */

  /** \defgroup linsol
  * @copydoc main_linsol
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_linsol
  * \endif
  */
  ///@{
  CASADI_EXPORT Function linsol(const std::string& name, const std::string& solver,
                                const Sparsity& sp, int nrhs, const Dict& opts=Dict());
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_linsol(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_linsol(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_linsol(const std::string& name);

  /** \brief Number of linear solver inputs */
  CASADI_EXPORT int linsol_n_in();

  /** \brief Number of linear solver outputs */
  CASADI_EXPORT int linsol_n_out();

  ///@{
  /** \brief Linear solver input scheme */
  CASADI_EXPORT std::vector<std::string> linsol_in();
  CASADI_EXPORT std::string linsol_in(int ind);
  ///@}

  ///@{
  /** \brief Linear solver output scheme */
  CASADI_EXPORT std::vector<std::string> linsol_out();
  CASADI_EXPORT std::string linsol_out(int ind);
  ///@}

  /** @} */

} // namespace casadi

#endif // CASADI_LINSOL_HPP
