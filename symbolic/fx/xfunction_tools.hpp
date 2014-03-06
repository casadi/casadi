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

#ifndef XFUNCTION_TOOLS_HPP
#define XFUNCTION_TOOLS_HPP

#include "fx.hpp"
#include "sx_function.hpp"
#include "mx_function.hpp"

namespace CasADi{

  /// @{
  /** \brief Make a vector-valued function out of a matrix-valued one.
   *  In spirit, this function is like applying vec() to all inputs outputs
   */
  SXFunction vec(SXFunction f);
  MXFunction vec(FX f);
  /// @}  
                
} // namespace CasADi


#endif // XFUNCTION_TOOLS_HPP
