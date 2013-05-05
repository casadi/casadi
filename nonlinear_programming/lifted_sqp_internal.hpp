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

#ifndef LIFTED_SQP_INTERNAL_HPP
#define LIFTED_SQP_INTERNAL_HPP

#include "lifted_sqp.hpp"
#include "symbolic/fx/nlp_solver_internal.hpp"
#include "symbolic/fx/qp_solver.hpp"

namespace CasADi{
    
class LiftedSQPInternal : public NLPSolverInternal{

public:
  explicit LiftedSQPInternal(const FX& F, const FX& G);
  virtual ~LiftedSQPInternal();
  virtual LiftedSQPInternal* clone() const{ return new LiftedSQPInternal(*this);}
  
  virtual void init();
  virtual void evaluate(int nfdir, int nadir);
  
  /// QP solver for the subproblems
  QPSolver qp_solver_;

  /// maximum number of sqp iterations
  int maxiter_; 

  /// stopping criterion for the stepsize
  double toldx_;
  
  /// stopping criterion for the lagrangian gradient
  double tolgl_;

  /// Linesearch parameters
  //@{
  double sigma_;
  double rho_;
  double mu_safety_;
  double eta_;
  double tau_;
  int maxiter_ls_;
  //@}
  
   /// Residual
   DMatrix d_k_;
   
   /// Primal step
   DMatrix dx_k_;
   
   /// Dual step
   DMatrix dlam_x_k_, dlam_g_k_;

   /// Indices
   enum GIn{G_X,G_LAM_X,G_LAM_G,G_NUM_IN};
   enum GOut{G_D,G_G,G_F,G_NUM_OUT};
   
   enum LinIn{LIN_X,LIN_LAM_X,LIN_LAM_G,LIN_D,LIN_NUM_IN};
   enum LinOut{LIN_F1,LIN_J1,LIN_F2,LIN_J2,LIN_NUM_OUT};

   enum ExpIn{ERR_X,EXP_EXP_X,EXP_LAM_G,EXP_D,EXP_DU,EXP_DLAM_F2,EXP_NUM_IN};
   enum ExpOut{EXP_E,EXP_NUM_OUT};

   /// Residual function
   FX rfcn_;
   
   /// Quadratic approximation
   FX lfcn_;
   
   /// Step expansion
   FX efcn_;
   
   /// Dimensions
   int nu, nv;
   
};

} // namespace CasADi

#endif //LIFTED_SQP_INTERNAL_HPP
