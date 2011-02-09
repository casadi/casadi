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

#ifndef CASADI_OCP_HPP
#define CASADI_OCP_HPP

#include <vector>

#include "fx.hpp"

namespace CasADi{

// Forward declaration of internal class
class OCPInternal;

/** \brief OCP execution of functions
  \author Joel Andersson
  \date 2011
*/ 
class OCP : public FX{
public:

  /// Default constructor
  OCP();

  /// Create a OCP
  explicit OCP(const std::vector<FX>& L, // Cost functions
               const std::vector<FX>& F, // Dynamic equation
               const std::vector<FX>& H=std::vector<FX>(), // Path constraints
               const std::vector<FX>& G=std::vector<FX>()); // Coupling constraints

  /// Access functions of the node
  OCPInternal* operator->();

  /// Const access functions of the node
  const OCPInternal* operator->() const;
  
};

} // namespace CasADi


#endif // CASADI_OCP_HPP

