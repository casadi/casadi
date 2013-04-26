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

#ifndef DERIVATIVE_HPP
#define DERIVATIVE_HPP

#include <vector>

#include "fx.hpp"

namespace CasADi{

// Forward declaration of internal class
class DerivativeInternal;

/** \brief Derivative class
        
  Universal directional derivative class, calculates a number of forward and adjoint
  directional derivatives of a function based on operator overloading AD.
  
  This is an internal class. Users should use the syntax f.derivative()
  
  \author Joel Andersson 
  \date 2012
*/ 
class Derivative : public FX{
  friend class FXInternal;
public:
  
  /// Default constructor
  Derivative();

  /// Create a Derivative
  explicit Derivative(const FX& fcn, int nfwd, int nadj);
  
  /// Access functions of the node
  DerivativeInternal* operator->();

  /// Const access functions of the node
  const DerivativeInternal* operator->() const;
  
};

} // namespace CasADi


#endif // DERIVATIVE_HPP

