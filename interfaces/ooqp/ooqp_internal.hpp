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
  *  Because constraints after initialization, we cannot allocate a Solver during init
  */
  virtual void allocate();
  
  virtual void evaluate(int nfdir, int nadir);
  
  protected:
    /** 
     * OOQP specific pointers.
     * @ {
     * */
    QpGenSparseMa27 * qp_;
    QpGenData      * prob_;
    QpGenVars      * vars_; 
    QpGenResiduals * resid_;
    GondzioSolver  * s_;
    /* @} */
    
    /// Boolean vector indicating equality with 0 and inequality with 1
    std::vector<int> constraints_;
    
    /// The non-zero indices of the equality constraints (vector of size n_eq)
    std::vector<int> eq_;

    /// The non-zero indices of the inequality constraints (vector of size n_eq)
    std::vector<int> ineq_;
    
    /// Number of equality constraints
    int n_eq_;
    
    /// Number of inequality constraints
    int n_ineq_;
    
    /// NZ indices of H elements  on the lower triangular side
    std::vector<int> nz_;
    
    /// The lower triangular part of H
    DMatrix Hl_;
    
    /// The non-zero indices into H that make up Hl 
    std::vector<int> Hl_nz_;
    
    /// The equality part of A
    DMatrix A_eq_;
    
    /// The inequality part of A
    DMatrix A_ineq_;
    
    /// The equality part of LBA
    DMatrix BA_eq_;
    
    /// The inequality part of LBA
    DMatrix LBA_ineq_;
    
    /// The inequality part of UBA
    DMatrix UBA_ineq_;
    
    /// The LBX with -inf substituted by 0 (needed for OOQP)
    DMatrix LBX_;
    
    /// The UBX with inf substituted by 0 (needed for OOQP)
    DMatrix UBX_;
    
    /// A vector to do slicing operations on A. 
    std::vector<int> all_A_;
    
    /** 
     *  In OOQP, infinite bounds need a special treatment. They must be deactivated by setting the i-something std::vector<char>.
     *  they have to be zero for infinite bounds and 1 otherwise
     * @{
     */
    std::vector<char> ixlow_;
    std::vector<char> ixupp_;

    std::vector<char> iclow_;
    std::vector<char> icupp_;
    
    /**
     * Because OOQP does not do const correctness properly, and because we cannot trust it, we copy all(rowind,col) data
     * @{
     * */
    std::vector<int> Hl_rowind_;
    std::vector<int> Hl_col_;
    
    std::vector<int> eq_rowind_;
    std::vector<int> eq_col_;
    
    std::vector<int> ineq_rowind_;
    std::vector<int> ineq_col_;
    /** @} */
    
    
    /// Throw error
    static void ooqp_error(const std::string& module, int flag);
    
    /// Calculate the error message map
    static std::map<int,std::string> calc_flagmap();
    
    /// Error message map
    static std::map<int,std::string> flagmap;
};


} // namespace Interfaces
} // namespace CasADi

#endif //OOQP_INTERNAL_HPP

