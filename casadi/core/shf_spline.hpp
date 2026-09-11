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

#ifndef CASADI_SHF_SPLINE_HPP
#define CASADI_SHF_SPLINE_HPP

#include "function.hpp"

namespace casadi {
  /** \brief Smooth Hat Function (SHF) spline with fixed table values
   *
   * Approximates the piecewise-linear interpolant of a look-up table by a
   * k-times continuously differentiable spline that coincides with the table
   * exactly outside epsilon-neighbourhoods of the interior grid points, and
   * converges to it as epsilon goes to zero. The coefficients are the table
   * values themselves; there is no fitting step.
   *
   * The resulting Function has inputs (x, eps) and one output of size m-by-1.
   * eps holds one absolute half-width per axis and must stay below half the
   * smallest grid spacing of the axis it applies to. Only x is differentiable.
   */
  CASADI_EXPORT Function shf_spline(const std::string& name,
    const std::vector< std::vector<double> >& grid,
    const std::vector<double>& values,
    casadi_int k, casadi_int m,
    const Dict& opts=Dict());

  /** \brief Smooth Hat Function (SHF) spline with parametric table values
   *
   * Like shf_spline with fixed values, but the table values arrive as an
   * additional (non-differentiable) input: inputs are (x, C, eps).
   */
  CASADI_EXPORT Function shf_spline(const std::string& name,
    const std::vector< std::vector<double> >& grid,
    casadi_int k, casadi_int m,
    const Dict& opts=Dict());

} // namespace casadi

#endif // CASADI_SHF_SPLINE_HPP
