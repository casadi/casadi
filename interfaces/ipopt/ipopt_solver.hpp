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

#ifndef IPOPT_SOLVER_HPP
#define IPOPT_SOLVER_HPP

#include "symbolic/fx/nlp_solver.hpp"

namespace CasADi{
  
class IpoptInternal;
  
// List from ipopt_internal.cpp
/**
* \brief interface to IPOPT NLP solver
* @copydoc NLPSolver_doc
*
* When in warmstart mode, output NLP_SOLVER_LAM_X may be used as input
*
* NOTE: Even when max_iter == 0,  it is not guaranteed that input(NLP_SOLVER_X0) == output(NLP_SOLVER_X). Indeed if bounds on X or constraints are unmet, they will differ.
*       
*  A good resource about the algorithms in IPOPT is: Wachter and L. T. Biegler, On the Implementation of an Interior-Point Filter Line-Search Algorithm for Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57, 2006 (As Research Report RC 23149, IBM T. J. Watson Research Center, Yorktown, USA
*
* Caveats: 
*   * with default options, multipliers for the decision variables are wrong for equality constraints.
*     Change the 'fixed_variable_treatment' to 'make_constraint' or 'relax_bounds' to obtain correct results.
* 
*
*/
class IpoptSolver : public NLPSolver {
  public:
    /// Default constructor
    IpoptSolver();

    /// \brief Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit IpoptSolver(const FX& F,         /**< F objective function: \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}]\f$*/
                         const FX& G = FX(),  /**< constraint function (default only bound constraints): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^m]\f$ */
                         const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory): \f$ [\mathbf{R}^n, \mathbf{R}^m, \mathbf{R}] \mapsto [\mathbf{R}^{n x n}]\f$ \n The third input argument for H is \f$ \sigma \f$, a scaling factor for F. */
                         const FX& J = FX(),  /**< Jacobian of G (default -> differentiate): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^{m x n}]\f$ */
                         const FX& GF = FX()  /**< Gradient of the objective function (default: adjoint mode AD on F): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^n]\f$ */
                        );

    FX getGF()	const;
    
    /** \brief Get the reduced Hessian. 
     * Requires a patched sIPOPT installation, see CasADi documentation. */
    DMatrix getReducedHessian();
    
    /// Access functions of the node
    IpoptInternal* operator->();
    const IpoptInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function 
    #ifdef SWIG
    %callback("%s_cb");
    #endif
    static NLPSolver creator(const FX& F, const FX& G, const FX& H, const FX& J){ return IpoptSolver(F,G,H,J);}
    #ifdef SWIG
    %nocallback;
    #endif

};

} // namespace CasADi

#endif //IPOPT_SOLVER_HPP
