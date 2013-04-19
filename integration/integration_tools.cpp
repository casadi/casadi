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

#include "integration_tools.hpp"

using namespace std;

namespace CasADi{
  
  std::vector<double> collocationPoints(int order, const std::string& scheme) {
    if (scheme=="radau") {
      casadi_assert_message(order>0 && order<10,"Error in collocationPoints(order, scheme): only order up to 9 supported for scheme 'radau', but got " << order << ".");
      return std::vector<double>(radau_points[order],radau_points[order]+order+1);
    } else if (scheme=="legendre") {
      casadi_assert_message(order>0 && order<10,"Error in collocationPoints(order, scheme): only order up to 9 supported for scheme 'legendre', but got " << order << ".");
      return std::vector<double>(legendre_points[order],legendre_points[order]+order+1);
    } else {
      casadi_error("Error in collocationPoints(order, scheme): unknown scheme '" << scheme << "'. Select one of 'radau', 'legendre'.");
    }
  }

} // namespace CasADi

