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

#ifndef KINSOL_INTERNAL_HPP
#define KINSOL_INTERNAL_HPP

#include "kinsol_solver.hpp"
#include "casadi/fx/implicit_function_internal.hpp"
#include "casadi/fx/linear_solver.hpp"
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_dense.h>  /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>  /* definition of type double */
#include <kinsol/kinsol.h>            /* prototypes for CVode fcts. and consts. */
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h> 
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_impl.h> /* Needed for the provided linear solver */
#include <ctime>

namespace CasADi{
namespace Sundials{
  
class KinsolInternal : public ImplicitFunctionInternal{
  friend class KinsolSolver;
public:
  /** \brief  Constructor */
  explicit KinsolInternal(const FX& f, int nrhs);

  /** \brief  Clone */
  virtual KinsolInternal* clone() const;

  /** \brief  Destructor */
  virtual ~KinsolInternal();

  /** \brief  Initialize stage */
  virtual void init();
  
  /** \brief  Evaluate */
  virtual void evaluate(int fsens_order, int asens_order);

  /** \brief Generate a linear solver for the sensitivity equations */
  KinsolSolver jac(int iind=0, int oind=0);

  /** \brief Generate a linear solver for the sensitivity equations */
  KinsolSolver jac(const std::vector<int> iind, int oind=0);

  /** \brief Residual */
  void func(N_Vector u, N_Vector fval);
  void djac(int N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2);

  /** \brief Wrappers */
  static int func_wrapper(N_Vector u, N_Vector fval, void *user_data);
  static int djac_wrapper(int N, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2);
  
  // KINSOL memory block
  void* mem_;
  
  // Jacobian
  FX J_;
  
  // Variable
  N_Vector u_;
  
  // Scaling
  N_Vector u_scale_, f_scale_;
  
  // For timings
  clock_t time1_, time2_;
  
  // Accummulated time since last reset:
  double t_func_; // time spent in the residual function
  double t_jac_; // time spent in the jacobian function

  // Globalization strategy
  int strategy_;

  // Nonlinear solver for the augmented system
  KinsolSolver aug_;
};


} // namespace Sundials
} // namespace CasADi

#endif //KINSOL_INTERNAL_HPP

