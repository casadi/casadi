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

#ifndef IP_METHOD_HPP
#define IP_METHOD_HPP

#include "casadi/fx/nlp_solver.hpp"

namespace CasADi{
  
class IPInternal;
  
/**
  \brief Interior point method
  This method is experimental only. Do not attempt to use if you do not intend to dive into the source code.
  The current purpose of the class is to show how an IP method can be implemeted in CasADi.
  If someone wants to take responsibility for this class and make it work, then please contact the CasADi developers.
  
  \author Joel Andersson
  \date 2012
*/
class IPMethod : public NLPSolver {
  public:
    /// Default constructor
    IPMethod();

    /// \brief Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit IPMethod(const FX& F,         /**< F objective function: \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}]\f$*/
                      const FX& G = FX(),  /**< constraint function (default only bound constraints): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^m]\f$ */
                      const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory): \f$ [\mathbf{R}^n, \mathbf{R}^m, \mathbf{R}] \mapsto [\mathbf{R}^{n x n}]\f$ \n The third input argument for H is \f$ \sigma \f$, a scaling factor for F. */
                      const FX& J = FX()   /**< Jacobian of G (default -> differentiate): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^{m x n}]\f$ */
                      );

    /// Access functions of the node
    IPInternal* operator->();
    const IPInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function 
    #ifdef SWIG
    %callback("%s_cb");
    #endif
    static NLPSolver creator(const FX& F, const FX& G, const FX& H, const FX& J){ return IPMethod(F,G,H,J);}
    #ifdef SWIG
    %nocallback;
    #endif

    /// @Joris: This would be an alternative
    static NLPSolverCreator getCreator(){return creator;}

    
};

} // namespace CasADi

#endif //IP_METHOD_HPP
