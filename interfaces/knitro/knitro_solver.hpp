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

#ifndef KNITRO_SOLVER_HPP
#define KNITRO_SOLVER_HPP

#include "casadi/fx/nlp_solver.hpp"

namespace CasADi{
  
class KnitroInternal;
  
class KnitroSolver : public NLPSolver {
  public:
    /// Default constructor
    KnitroSolver();
    
    /// Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit KnitroSolver(const FX& F,         /**< F objective function */
                         const FX& G = FX(),  /**< constraint function (default only bound constraints) */
                         const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory) */
                         const FX& J = FX(),  /**< Jacobian of G (default -> differentiate) */
                         const FX& GF = FX()  /**< Gradient of the objective function (default: adjoint mode AD on F) */
                        );

    /// Access functions of the node
    KnitroInternal* operator->();
    const KnitroInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    
};

} // namespace CasADi

#endif //KNITRO_SOLVER_HPP
