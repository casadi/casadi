/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_CUSTOM_FUNCTION_INTERNAL_HPP
#define CASADI_CUSTOM_FUNCTION_INTERNAL_HPP

#include "custom_function.hpp"
#include "function_internal.hpp"
#include "../functor.hpp"
#include <string>

/// \cond INTERNAL

namespace casadi {

  /** \brief  Internal class for CustomFunction
  \author Joel Andersson
  \date 2010
  A regular user should never work with any Node class. Use CustomFunction directly.
  */
class CASADI_CORE_EXPORT CustomFunctionInternal : public FunctionInternal {
  friend class CustomFunction;
  public:

    /** \brief  Create a function */
    explicit CustomFunctionInternal(const CustomEvaluate &c_fcn,
                                    const std::vector<casadi::Sparsity> &inputscheme,
                                    const std::vector<casadi::Sparsity> &outputscheme);

    /** \brief  Destructor */
    virtual ~CustomFunctionInternal();

    /** \brief  Cloning */
    virtual CustomFunctionInternal* clone() const { return new CustomFunctionInternal(*this);}

    /** \brief  Evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    CustomEvaluate evaluate_;

    /// A reference to this object to be passed to the user functions
    CustomFunction ref_;

}; // class CustomFunctionInternal



} // namespace casadi
/// \endcond

#endif // CASADI_CUSTOM_FUNCTION_INTERNAL_HPP
