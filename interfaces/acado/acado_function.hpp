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

#ifndef ACADO_FUNCTION_HPP
#define ACADO_FUNCTION_HPP

#include <casadi/fx/fx.hpp>
#include "acado_forward_declarations.hpp"

namespace CasADi{

/** \brief  CasADi to ACADO function interface */
class AcadoFunction{
  public:
      // Constructor
      AcadoFunction(const FX& f=FX());
      
      // Destructor
      ~AcadoFunction();
      
      // Initialize
      void init();

      // CasADi function
      FX f_;
      
      // ACADO c function pointer
      ACADO::CFunction *fcn_;

      // Input dimensions
      std::vector<int> dim_;

      // Number of equations
      int neq_;
      
      // C-callback functions called by ACADO
      static void fcn_wrapper( double *x, double *res, void *user_data );
      static void fcn_fwd_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data);
      static void fcn_adj_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data);

      // Internal functions
      void fcn( double *x, double *res);
      void fcn_fwd(int number, double *x, double *seed, double *f, double *df);
      void fcn_adj(int number, double *x, double *seed, double *f, double *df);
};

} // namespace CasADi

#endif //ACADO_FUNCTION_HPP
