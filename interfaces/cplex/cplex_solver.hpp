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
#ifndef CPLEX_SOLVER_HPP
#define CPLEX_SOLVER_HPP

#include "casadi/fx/nlp_solver.hpp"

namespace CasADi{
  
class CplexInternal;
  
/** \brief Interface to CPLEX solver.
  @copydoc NLPSolver_doc
  Attention! The interface is not complete yet.
  Also if a quadratic term can be set with this interface, it is ignored!
  \author Carlo Savorgnan
  \date 2011
*/
class CplexSolver : public NLPSolver {
  // TODO comment me!!!!
  public:
    /// Default constructor
    CplexSolver();
    
    /// Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit CplexSolver(const FX& F,         /**< F objective function */
                         const FX& G = FX(),  /**< constraint function (default only bound constraints) */
                         const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory). NOT USED*/
                         const FX& J = FX(),  /**< Jacobian of G (default -> differentiate) */
                         const FX& GF = FX()  /**< Gradient of the objective function (default: adjoint mode AD on F) */
                        );

    /// Access functions of the node
    CplexInternal* operator->();
    const CplexInternal* operator->() const;
    
    /// Set CPLEX integer parameters
    void setIntParam(const std::string& name, int val);
    
    /// Set CPLEX double parameters
    void setDoubleParam(const std::string& name, double val);

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    
};

} // namespace CasADi

#endif //CPLEX_SOLVER_HPP
