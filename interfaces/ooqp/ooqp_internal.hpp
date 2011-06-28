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

#ifndef OOQP_INTERNAL_HPP
#define OOQP_INTERNAL_HPP

#include "casadi/fx/qp_solver_internal.hpp"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "GondzioSolver.h"
#include "QpGenSparseMa27.h"

namespace CasADi{
namespace Interfaces {
  
class OOQPInternal : public QPSolverInternal {
  friend class OOQPSolver;
public:
  /** \brief  Constructor */
  explicit OOQPInternal();

  /** \brief  Clone */
  virtual OOQPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit OOQPInternal(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A);

  /** \brief  Destructor */
  virtual ~OOQPInternal();

  /** \brief  Initialize */
  virtual void init();
  
  /** \brief  Allocate
  *
  *  Because constraints after initialization, we cannot allocate a Solver during init
  */
  virtual void allocate();
  
  virtual void evaluate(int nfdir, int nadir);
  
  protected:
    QpGenSparseMa27 * qp;
    QpGenData      * prob;
    QpGenVars      * vars; 
    QpGenResiduals * resid;
    GondzioSolver  * s;
    
    /// Boolean vector indicating equality with 0 and inequality with 1
    std::vector<int> constraints;
    
    /// The non-zero indices of the equality constraints (vector of size n_eq)
    std::vector<int> eq;

    /// The non-zero indices of the inequality constraints (vector of size n_eq)
    std::vector<int> ineq;
    
    /// Number of equality constraints
    int n_eq;
    
    /// Number of inequality constraints
    int n_ineq;
    
    /** \brief Sorts the constraints
    * \pre
    *   constraints member is set
    * \post 
    *   eq, ineq, n_eq, n_ineq  are set
    *   ineq.size()==n_ineq,  eq.size()==n_eq
    */
    void sort_constraints();
    
    /** \brief Guess the constraints from QP_LBA and QP_UBA if 'constraint' option is not set
    * \post 
    *   constraints member is set
    */
    void guess_constraints();
    
    /// NZ indices of H elements  on the lower triangular side
    std::vector<int> nz;
    
    /// The lower triangular part of H
    DMatrix Hl;
    
    /// The equality part of A
    DMatrix A_eq;
    
    /// The inequality part of A
    DMatrix A_ineq;
    
    /// The equality part of LBA
    DMatrix BA_eq;
    
    /// The inequality part of LBA
    DMatrix LBA_ineq;
    
    /// The inequality part of UBA
    DMatrix UBA_ineq;
    
    //xlow, ixlow are the lower bounds on x. These contain the information in the
    //lower bounding vector l in the formulation given above. 
    //If there is a bound on element k of x (that is, lk > -1), then xlow[k] should
    //be set to the value of lk and ixlow[k] should be set to one. 
    std::vector<int> ixlow;
    std::vector<int> xlow;

    std::vector<int> ixupp;
    std::vector<int> xupp;
    
    std::vector<int> iclow;
    std::vector<int> icupp;
    
};


} // namespace Interfaces
} // namespace CasADi

#endif //OOQP_INTERNAL_HPP

