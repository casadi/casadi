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


#ifndef CASADI_KERNEL_SUM_2D_HPP
#define CASADI_KERNEL_SUM_2D_HPP

#include "function.hpp"

namespace casadi {

  /** \brief  Forward declaration of internal class */
  class KernelSum2DInternal;

  /** KernelSum2D 
  
        Consider a dense matrix V.
        
        KernelSum computes     
  
        F(V,X)  = sum_i sum_j  f ( [i;j], V(i,j), X)
        
          with X: [x;y]
          
        where the summation is taken for all entries (i,j)
        that are a distance r away from X.
        
        This function assumes that V is fixed: 
        sensitivities with respect to it are not computed.
        
        This allows for improved speed of evaluation.
        
        Having V fixed is a common use case:
          V may be a large bitmap (observation),
          onto which a kernel is fitted.
          
        Summation does not occur outside the image.
        Runtime will not grow after distance r grows large enough to contian the whole image.
        
  
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT KernelSum2D : public Function {
  public:
    /** \brief Default constructor */
    KernelSum2D();

    /** \brief Constructor (generic kernel_sum_2d) */
    KernelSum2D(const std::string& name, const Function& f,
           const std::pair<int, int> & size,
           double r,
           int n,
           const Dict& opts=Dict());


    /** \brief  Access functions of the node */
    KernelSum2DInternal* operator->();

    /** \brief  Const access functions of the node */
    const KernelSum2DInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };

} // namespace casadi

#endif // CASADI_KERNEL_SUM_2D_HPP
