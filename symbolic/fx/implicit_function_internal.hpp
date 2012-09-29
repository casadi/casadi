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

#ifndef IMPLICIT_FUNCTION_INTERNAL_HPP
#define IMPLICIT_FUNCTION_INTERNAL_HPP

#include "implicit_function.hpp"
#include "fx_internal.hpp"

namespace CasADi{
  
// Forward declaration of internal class
class ImplicitFunctionInternal;

/// Internal class
class ImplicitFunctionInternal : public FXInternal{
  public:
    /** \brief Constructor
    *
    * \param f   FX mapping from (n+1) inputs to 1 output.
    */
    ImplicitFunctionInternal(const FX& f, int nrhs);
        
    /// Destructor
    virtual ~ImplicitFunctionInternal() = 0;
    
    /// Initialize
    virtual void init();
    
    /// Solve the system of equations
    virtual void evaluate(int nfdir, int nadir) = 0;
    
    /// Sparsity in CRS format
    FX f_;
    
    /// Number of equations
    int N_;
    
    /// Number of right hand sides
    int nrhs_;
    
    /// Number of forward derivative directions of the function
    int nfdir_fcn_;
    
    /// Number of adjoint derivative directions of the function
    int nadir_fcn_;
};



} // namespace CasADi

#endif //IMPLICIT_FUNCTION_INTERNAL_HPP

