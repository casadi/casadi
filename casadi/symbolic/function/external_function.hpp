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

#ifndef EXTERNAL_FUNCTION_HPP
#define EXTERNAL_FUNCTION_HPP

#include "function.hpp"
#include <string>

namespace casadi{


/** \brief  Forward declaration of internal class */
class ExternalFunctionInternal;

/** \brief  Interface for a function that is not implemented by CasADi symbolics
  \author Joel Andersson
  \date 2011
  */
class CASADI_SYMBOLIC_EXPORT ExternalFunction : public Function{

public:

/** \brief  CONSTRUCTORS: */
/** \brief  default constructor */
  ExternalFunction();

  /** \brief  Create an empty function */
  explicit ExternalFunction(const std::string& bin_name);

  /** \brief  Access functions of the node */
  ExternalFunctionInternal* operator->();

  /** \brief  Const access functions of the node */
  const ExternalFunctionInternal* operator->() const;

  /** \brief  Check if the pointer points towards a valid object */
  virtual bool checkNode() const;

};

} // namespace casadi


#endif // EXTERNAL_FUNCTION_HPP
