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


#ifndef CASADI_ORACLE_FUNCTION_HPP
#define CASADI_ORACLE_FUNCTION_HPP

#include "function_internal.hpp"
#include "../timing.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Function memory with temporary work vectors */
  struct CASADI_EXPORT OracleMemory {
    // Work vectors
    const double** arg;
    double** res;
    int* iw;
    double* w;

    // Function specific statistics
    std::map<std::string, FStats> fstats;
  };

  /** \brief Base class for functions that perform calculation with an oracle
      \author Joel Andersson
      \date 2016
  */
  class CASADI_EXPORT OracleFunction : public FunctionInternal {
  public:
    /** \brief  Constructor */
    OracleFunction(const std::string& name, const Function& oracle);

    /** \brief  Destructor */
    virtual ~OracleFunction() = 0;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_ORACLE_FUNCTION_HPP
