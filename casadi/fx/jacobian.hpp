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

#ifndef JACOBIAN_HPP
#define JACOBIAN_HPP

#include <vector>

#include "fx.hpp"

namespace CasADi{

// Forward declaration of internal class
class JacobianInternal;

/** \brief Jacobian class
	
  Universal Jacobian class, calculates the Jacobian of a function based on AD forward or adjoint.
  
  Options:\n
  "finite_differences" false\n
  "ad_mode"            "forward", "adjoint" or "default", i.e. forward if n_<=m_, otherwise adjoint\n
  "sparse"             false\n
  
  Any CasADi::FX can be used to take the Jacobian of.
  
  If the iind'th input argument has shape (m,n) and the oind'th output argument has shape (k,l),
  the output of this Jacobian will be of shape (k*l) x (m*n) .
  
  \author Joel Andersson 
  \date 2010
*/ 
class Jacobian : public FX{
public:

  /// Default constructor
  Jacobian();

  /// Create a Jacobian
  explicit Jacobian(const FX& fcn, int iind=0, int oind=0);

  /// Create a set of Jacobians (new formulation)
  Jacobian(const FX& fcn, const std::vector<std::pair<int,int> >& jblocks);
  
  /// Access functions of the node
  JacobianInternal* operator->();

  /// Const access functions of the node
  const JacobianInternal* operator->() const;
  
};

} // namespace CasADi


#endif // JACOBIAN_HPP

