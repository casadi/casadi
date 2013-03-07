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
#include "symbolic/fx/linear_solver.hpp"

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
    
    /** \brief  Update the number of sensitivity directions during or after initialization */
    virtual void updateNumSens(bool recursive);

    /// Solve the system of equations
    virtual void evaluate(int nfdir, int nadir) = 0;
    
    /// The function F(z, x1, x2, ..., xn) == 0
    FX f_;
    
    /// Number of equations
    int N_;
    
    /// Number of right hand sides
    int nrhs_;
    
    /** \brief  Create a new ImplicitFunctionInternal */
    virtual ImplicitFunctionInternal* create(const FX& f, int nrhs=1) const = 0;
    
    /** \brief Generate a linear solver for the sensitivity equations */
    ImplicitFunction jac(int iind, int oind=0);
    
    /** \brief Generate a linear solver for the sensitivity equations */
    ImplicitFunction jac(const std::vector<int> iind, int oind=0);
    
    /// Set the jacobian of F
    void setJacobian(FX &J);
    
    /// Jacobian
    FX J_;
    
    /// Linear solver
    LinearSolver linsol_; 
  protected:
  
    /** Calculate sensitivities of implicit solver
    * \param linsol_prepared may specify that the linear solver is already prepared with the Jacobian
    */
    void evaluate_sens(int nfdir, int nadir, bool linsol_prepared=false);
};



} // namespace CasADi

#endif //IMPLICIT_FUNCTION_INTERNAL_HPP

