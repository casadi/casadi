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

#ifndef C_FUNCTION_HPP
#define C_FUNCTION_HPP

#include "fx.hpp"
#include <string>

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class CFunctionInternal;

// Forward declaration
class CFunction;

//#ifndef SWIG
/** \brief  Wrapper around functions
 *
 * For wrappers around Python functions, see the PyFunction example at the bottom.
 */
typedef void (*CFunctionWrapper)(CFunction &f, int nfdir, int nadir, void* user_data);
//#endif // SWIG

/** \brief  Interface to function implemented as plain code 
  \author Joel Andersson 
  \date 2010
*/
class CFunction : public FX{

public:

/** \brief  default constructor */
  CFunction();

  //@{
  /** \brief  Create a function with input/output schemes given */
  explicit CFunction(CFunctionWrapper c_fcn, const std::vector<CRSSparsity> &inputscheme, const std::vector<CRSSparsity> &outputscheme);
  
  explicit CFunction(CFunctionWrapper c_fcn, const IOSchemeVector< CRSSparsity > &inputscheme, const std::vector<CRSSparsity> &outputscheme);
  
  explicit CFunction(CFunctionWrapper c_fcn, const std::vector<CRSSparsity> &inputscheme, const IOSchemeVector< CRSSparsity > &outputscheme);
  
  explicit CFunction(CFunctionWrapper c_fcn, const IOSchemeVector< CRSSparsity > &inputscheme, const IOSchemeVector< CRSSparsity > &outputscheme);
  //@}
  
  /** \brief  Create a function, user sets inputs outputs manually */
  explicit CFunction(CFunctionWrapper c_fcn);

  /** \brief  Access functions of the node */
  CFunctionInternal* operator->();
  
  /** \brief  Const access functions of the node */
  const CFunctionInternal* operator->() const;
  
  /** \brief  Check if the pointer points towards a valid object */
  virtual bool checkNode() const;

}; // class CFunction
  


} // namespace CasADi


#endif // C_FUNCTION_HPP
