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


#ifndef CASADI_EXPM_HPP
#define CASADI_EXPM_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_expm Title
      \par

      Performs a matrix exponentiation
      expm(A)


      \generalsection{Expm}
      \pluginssection{Expm}

      \author Joris Gillis
      \date 2017

      \identifier{21l} */

  /** \defgroup expm Title
  * @copydoc main_expm
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_expm
  * \endif
  */
  ///@{
  CASADI_EXPORT Function expmsol(const std::string& name, const std::string& solver,
                           const Sparsity& A, const Dict& opts=Dict());
  ///@}

  /** \brief Get the number of expm solver inputs

      \identifier{rs} */
  CASADI_EXPORT casadi_int expm_n_in();

  /** \brief Get the number of expm solver outputs

      \identifier{rt} */
  CASADI_EXPORT casadi_int expm_n_out();

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_expm(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_expm(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_expm(const std::string& name);

  /** @} */
} // namespace casadi

#endif // CASADI_EXPM_HPP
