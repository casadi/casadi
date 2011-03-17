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

#ifndef X_FUNCTION_HPP
#define X_FUNCTION_HPP

#include "fx.hpp"
#include "../sx/sx.hpp"
#include <vector>

namespace CasADi{

/// Forward declaration of internal class
class XFunctionInternal;

/**   \brief Dynamically created function that can be expanded into a series of scalar operations. Base class for XFunction and MXFunction.
\author Joel Andersson 
\date 2011
*/

class XFunction : public FX{

public:
  /// Default constructor
  XFunction();

  /// Access functions of the node 
  XFunctionInternal* operator->();

  /// Const access functions of the node 
  const XFunctionInternal* operator->() const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
};

#ifdef SWIG
%extend XFunction{
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif

} // namespace CasADi

#endif // X_FUNCTION_HPP
