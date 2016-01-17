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


#ifndef CASADI_ROOTFINDER_HPP
#define CASADI_ROOTFINDER_HPP

#include "function.hpp"
#include "linsol.hpp"

namespace casadi {

  /** \defgroup main_rootfinder
   * Create a solver for rootfinding problems
   * Takes a function where one of the inputs is unknown and one of the outputs
   * is a residual function that is always zero, defines a new function where
   * the the unknown input has been replaced by a _guess_ for the unknown and the
   * residual output has been replaced by the calculated value for the input.
   *
   * For a function
   * [y0, y1, ...,yi, .., yn] = F(x0, x1, ..., xj, ..., xm),
   * where xj is unknown and yi=0, defines a new function
   * [y0, y1, ...,xj, .., yn] = G(x0, x1, ..., xj_guess, ..., xm),
   *
   * xj and yi must have the same dimension and d(yi)/d(xj) must be invertable.
   *
   * By default, the first input is unknown and the first output is the residual.
   *
   *
   * \generalsection{Rootfinder}
   * \pluginssection{Rootfinder}
   *
   * \author Joel Andersson
   * \date 2011-2015
   */


  /** \defgroup rootfinder
  * @copydoc main_rootfinder
  *  @{
  */


  /** \if EXPANDED
  * @copydoc main_rootfinder
  * \endif
  */
  ///@{
  CASADI_EXPORT Function rootfinder(const std::string& name, const std::string& solver,
                               const Function& f, const Dict& opts=Dict());
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_rootfinder(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_rootfinder(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_rootfinder(const std::string& name);

  /** @} */

} // namespace casadi

#endif // CASADI_ROOTFINDER_HPP
