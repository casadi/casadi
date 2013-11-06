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

#ifndef CUSTOM_FUNCTION_HPP
#define CUSTOM_FUNCTION_HPP

#include "fx.hpp"
#include <string>

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class CustomFunctionInternal;

// Forward declaration
class CustomFunction;

/** \brief  Interface to a custom function
  
  Note: max_number_of_fwd_dir and max_number_of_adj_dir will be default zero
  
  
  \author Joel Andersson 
  \date 2010
*/
class CustomFunction : public FX{

public:

/** \brief  default constructor */
  CustomFunction();

  //@{
  /** \brief  Create a function with input/output schemes given */
  explicit CustomFunction(const CustomEvaluate &c_fcn, const std::vector<CRSSparsity> &inputscheme, const std::vector<CRSSparsity> &outputscheme);
  
  explicit CustomFunction(const CustomEvaluate &c_fcn, const IOSchemeVector< CRSSparsity > &inputscheme, const std::vector<CRSSparsity> &outputscheme);
  
  explicit CustomFunction(const CustomEvaluate &c_fcn, const std::vector<CRSSparsity> &inputscheme, const IOSchemeVector< CRSSparsity > &outputscheme);
  
  explicit CustomFunction(const CustomEvaluate &c_fcn, const IOSchemeVector< CRSSparsity > &inputscheme, const IOSchemeVector< CRSSparsity > &outputscheme);
  //@}
  
  /** \brief  Create a function, user sets inputs outputs manually */
  explicit CustomFunction(const CustomEvaluate &c_fcn);

  /** \brief  Access functions of the node */
  CustomFunctionInternal* operator->();
  
  /** \brief  Const access functions of the node */
  const CustomFunctionInternal* operator->() const;
  
  /** \brief  Check if the pointer points towards a valid object */
  virtual bool checkNode() const;

}; // class CustomFunction
  


} // namespace CasADi


#endif // CUSTOM_FUNCTION_HPP
