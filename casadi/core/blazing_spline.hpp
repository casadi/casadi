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


#ifndef CASADI_BLAZING_SPLINE_FUNCTION_HPP
#define CASADI_BLAZING_SPLINE_FUNCTION_HPP

#include "function.hpp"

namespace casadi {
  /** \brief Construct a specialized parametric BSpline
   * 
   * Specialized in these ways:
   *  - order 3 is assumed
   *  - up to dimension 3 supported
   *  - a single scalar output (m=1)
   *

      \identifier{2b9} */
  CASADI_EXPORT Function blazing_spline(const std::string& name,
    const std::vector< std::vector<double> >& knots,
    const Dict& opts=Dict());

} // namespace casadi

#endif // CASADI_BLAZING_SPLINE_FUNCTION_HPP
