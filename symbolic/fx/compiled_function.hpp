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

#ifndef COMPILED_FUNCTION_HPP
#define COMPILED_FUNCTION_HPP

#include "fx.hpp"
#include <string>

namespace CasADi{

  
/** \brief  Forward declaration of internal class */
class CompiledFunctionInternal;

/** \brief  Interface for a function that is not implemented by CasADi symbolics
  \author Joel Andersson 
  \date 2011
  */
class CompiledFunction : public FX{

public:

/** \brief  CONSTRUCTORS: */
/** \brief  default constructor */
  CompiledFunction();

  /** \brief  Create an empty function */
  explicit CompiledFunction(const std::string& bin_name);

  /** \brief  Access functions of the node */
  CompiledFunctionInternal* operator->();
  
  /** \brief  Const access functions of the node */
  const CompiledFunctionInternal* operator->() const;
    
  /** \brief  Check if the pointer points towards a valid object */
  virtual bool checkNode() const;

};

} // namespace CasADi


#endif // COMPILED_FUNCTION_HPP
