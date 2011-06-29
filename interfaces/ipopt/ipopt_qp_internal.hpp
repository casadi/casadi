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

#ifndef IPOPT_QP_INTERNAL_HPP
#define IPOPT_QP_INTERNAL_HPP

#include "casadi/fx/qp_solver_internal.hpp"
#include "interfaces/ipopt/ipopt_solver.hpp"

namespace CasADi{
namespace Interfaces {

  /** \brief Internal class for IpoptQPSolver
   * 
   * */
class IpoptQPInternal : public QPSolverInternal {
  friend class IpoptQPSolver;
public:
  /** \brief  Constructor */
  explicit IpoptQPInternal();

  /** \brief  Clone */
  virtual IpoptQPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit IpoptQPInternal(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A);

  /** \brief  Destructor */
  virtual ~IpoptQPInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate(int nfdir, int nadir);
  
  protected:
    IpoptSolver solver;
    
    // an MX that represents H, but dependant on a common MX.
    MX H_;
    // an MX that represents G, but dependant on a common MX.
    MX G_;
    // an MX that represents A, but dependant on a common MX.
    MX A_;
};


} // namespace Interfaces
} // namespace CasADi

#endif //IPOPT_QP_INTERNAL_HPP

