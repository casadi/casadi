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


#ifndef CASADI_BLAS_HPP
#define CASADI_BLAS_HPP

#include <casadi/core/casadi_export.h>
#include <string>

namespace casadi {

  /** \defgroup main_blas Title
      \par

      Pluggable dense matrix-multiply backend. The built-in "reference"
      plugin requires no external dependency; "classic" and "blasfeo"
      plugins (when built) load lazily from their respective shared
      libraries on first use.

      \generalsection{Blas}
      \pluginssection{Blas}

      \author Joris Gillis
      \date 2026

      \identifier{2gb} */

  /* \brief Check if a particular BLAS plugin is available */
  CASADI_EXPORT bool has_blas(const std::string& name);

  /* \brief Explicitly load a BLAS plugin dynamically */
  CASADI_EXPORT void load_blas(const std::string& name);

  /* \brief Get the documentation string for a BLAS plugin */
  CASADI_EXPORT std::string doc_blas(const std::string& name);

} // namespace casadi


#endif // CASADI_BLAS_HPP
