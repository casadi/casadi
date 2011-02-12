/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef VARIABLE_TOOLS_HPP
#define VARIABLE_TOOLS_HPP

#include "variable.hpp"

namespace CasADi{
namespace OptimalControl{
  /// Get a vector of expressions from a vector of variables
  std::vector<SX> sx(const std::vector<Variable> v);

  /// Get a vector of derivative expressions from a vector of variables
  std::vector<SX> der(const std::vector<Variable> v);

  /// Get a vector of the nominal values of a vector of variables
  std::vector<double> nominal(const std::vector<Variable> v);

} // namespace OptimalControl
} // namespace CasADi

#endif // VARIABLE_TOOLS_HPP

