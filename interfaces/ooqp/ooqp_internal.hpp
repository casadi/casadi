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
    
    /// NZ indices of H elements  on the lower triangular side
    std::vector<int> nz;
    
    /// The lower triangular part of H
    DMatrix Hl;
    
    /// The non-zero indices into H that make up Hl 
    std::vector<int> Hl_nz;
    
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
    
    
    /// The LBX with inf substituted
    DMatrix LBX;
    
    /// The UBX with inf substituted
    DMatrix UBX;
    
    std::vector<int> all_A;
    
    //xlow, ixlow are the lower bounds on x. These contain the information in the
    //lower bounding vector l in the formulation given above. 
    //If there is a bound on element k of x (that is, lk > -1), then xlow[k] should
    //be set to the value of lk and ixlow[k] should be set to one. 
    std::vector<char> ixlow;
    //std::vector<double> xlow;

    std::vector<char> ixupp;
    //std::vector<double> xupp;
    
    std::vector<char> iclow;
    std::vector<char> icupp;
    
    std::vector<int> Hl_rowind;
    std::vector<int> Hl_col;
    
    std::vector<int> eq_rowind;
    std::vector<int> eq_col;
    
    std::vector<int> ineq_rowind;
    std::vector<int> ineq_col;
    
    
    // Throw error
    static void ooqp_error(const std::string& module, int flag);
    
    // Calculate the error message map
    static std::map<int,std::string> calc_flagmap();
    
    // Error message map
    static std::map<int,std::string> flagmap;
};


} // namespace Interfaces
} // namespace CasADi

#endif //OOQP_INTERNAL_HPP

