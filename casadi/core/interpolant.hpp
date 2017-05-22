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


#ifndef CASADI_INTERPOLANT_HPP
#define CASADI_INTERPOLANT_HPP

#include "function.hpp"

namespace casadi {


  /** \defgroup main_interpolant
   * An interpolant function for lookup table data
   *
   *  \generalsection{Interpolant}
   *  \pluginssection{Interpolant}
   * \author Joel Andersson
   * \date 2016
   */

  /** \defgroup interpolant
  * @copydoc main_interpolant
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_interpolant
  * \endif
  */
  ///@{
  CASADI_EXPORT Function interpolant(const std::string& name,
                                     const std::string& solver,
                                     const std::vector<std::vector<double> >& grid,
                                     const std::vector<double>& values,
                                     const Dict& opts=Dict());
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_interpolant(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_interpolant(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_interpolant(const std::string& name);

  /** @} */

} // namespace casadi

#endif // CASADI_INTERPOLANT_HPP
