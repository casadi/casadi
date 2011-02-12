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

#ifndef CASADI_OCP2_HPP
#define CASADI_OCP2_HPP

#include <vector>

#include "fx.hpp"

namespace CasADi{

enum OCP2Input{
  OCP2_T,                      // Time grid
  OCP2_X, OCP2_LBX, OCP2_UBX,    // Differential state with bounds
  OCP2_Z, OCP2_LBZ, OCP2_UBZ,    // Algebraic state with bounds
  OCP2_XP, OCP2_LBXP, OCP2_UBXP, // State deriatives with bounds
  OCP2_U, OCP2_LBU, OCP2_UBU,    // Controls with bounds
  OCP2_P, OCP2_LBP, OCP2_UBP,    // Parameters with bounds
  OCP2_LBH, OCP2_UBH,           // Bounds for the point constraints
  OCP2_LBG, OCP2_UBG,           // Bounds for the coupling constraints
  OCP2_NUM_IN};
  
// Forward declaration of internal class
class OCP2Internal;

/** \brief Optimal control problem formulation
  \author Joel Andersson
  \date 2011
*/ 
class OCP2 : public FX{
public:

  /// Default constructor
  OCP2();

  /// Create a OCP2
  explicit OCP2(const std::vector<FX>& L, // Cost functions
               const std::vector<FX>& F, // Dynamic equation
               const std::vector<FX>& H=std::vector<FX>(), // Path constraints
               const std::vector<FX>& G=std::vector<FX>()); // Coupling constraints

  /// Access functions of the node
  OCP2Internal* operator->();

  /// Const access functions of the node
  const OCP2Internal* operator->() const;
  
};

} // namespace CasADi


#endif // CASADI_OCP2_HPP

