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

#ifndef SQP_METHOD_HPP
#define SQP_METHOD_HPP

#include "symbolic/fx/nlp_solver.hpp"

namespace CasADi{
  
class SQPInternal;
  
/**
  \brief Sequential Quadratic Programming method.
  The algorithm is a classical SQP method with either exact (may be also provided) or 
  damped BFGS Lagrange Hessian approximation.
  Two different line-search algorithms are available.
  First, Armijo (Wolfe) condition with backtracking (suffers from Maratos effect).
  Second, a line-search method that checks if the merit function is lower
  than the last k values (no Maratos effect).
  Both methods employ the L1 merit function.
  
  The method solves the problems of form:
  \verbatim
  min          F(x)
  x
  
  subject to
            LBG <= G(x) <= UBG
            LBX <=   x  <= UBX
  \endverbatim
  
  Nonlinear equalities can be introduced by setting LBG and UBG equal at the correct positions.
  
  The method is still under development and should be used with cared
  
  \author Joel Andersson and Attila Kozma
  \date 2012
*/
class SQPMethod : public NLPSolver {
  public:
    /// Default constructor
    SQPMethod();

    /// \brief Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit SQPMethod(const FX& F,         /**< F objective function: \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}]\f$*/
                       const FX& G = FX(),  /**< constraint function (default only bound constraints): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^m]\f$ */
                       const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory): \f$ [\mathbf{R}^n, \mathbf{R}^m, \mathbf{R}] \mapsto [\mathbf{R}^{n x n}]\f$ \n The third input argument for H is \f$ \sigma \f$, a scaling factor for F. */
                       const FX& J = FX()   /**< Jacobian of G (default -> differentiate): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^{m x n}]\f$ */
                       );

    /// Access functions of the node
    SQPInternal* operator->();
    const SQPInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function 
    #ifdef SWIG
    %callback("%s_cb");
    #endif
    static NLPSolver creator(const FX& F, const FX& G, const FX& H, const FX& J){ return SQPMethod(F,G,H,J);}
    #ifdef SWIG
    %nocallback;
    #endif

    /// Access the QPSolver used internally
    const QPSolver getQPSolver() const;
    
};

} // namespace CasADi

#endif //SQP_METHOD_HPP
