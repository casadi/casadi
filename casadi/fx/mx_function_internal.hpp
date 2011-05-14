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

#ifndef MX_FUNCTION_INTERNAL_HPP
#define MX_FUNCTION_INTERNAL_HPP

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "mx_function.hpp"
#include "x_function_internal.hpp"
#include "../mx/mx_node.hpp"

namespace CasADi{

/** \brief  Internal node class for MXFunction
  \author Joel Andersson 
  \date 2010
*/
class MXFunctionInternal : public XFunctionInternal{
  friend class MXFunction;
  
  public:

    /** \brief  Multiple input, multiple output constructor, only to be accessed from MXFunction, therefore protected */
    MXFunctionInternal(const std::vector<MX>& input, const std::vector<MX>& output);

    /** \brief  Make a deep copy */
    virtual MXFunctionInternal* clone() const;

    /** \brief  Destructor */
    virtual ~MXFunctionInternal();

    /** \brief  Evaluate the algorithm */
    virtual void evaluate(int nfdir, int nadir);

    /** \brief  Print description */
    virtual void print(std::ostream &stream) const;

    /** \brief  Initialize */
    virtual void init();
    
    /** \brief Set the lifting function */
    void setLiftingFunction(LiftingFunction liftfun, void* user_data);

    /** \brief Jacobian via source code transformation (identity matrix seed in a particular direction) */
    std::vector<MX> jac(int iind);
    
    /** \brief Jacobian via source code transformation */
    std::vector<MX> jac(const std::vector<MX>& fseed);

    /** \brief  An elemenent of the algorithm, namely an MX node */
    typedef MXAlgEl AlgEl;

    /** \brief  All the runtime elements in the order of evaluation */
    std::vector<AlgEl> alg;
  
    /** \brief  Matrix expressions that are to be evaluated */
    std::vector<MX> outputv;
    std::vector<int> outputv_ind;

    /** \brief  Dependent expressions */
    std::vector<MX> inputv;
    std::vector<int> inputv_ind;
      
    // Lifting function
    LiftingFunction liftfun_;
    void* liftfun_ud_;
};

} // namespace CasADi


#endif // MX_FUNCTION_INTERNAL_HPP

