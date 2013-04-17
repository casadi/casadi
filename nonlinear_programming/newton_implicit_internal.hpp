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

#ifndef NEWTON_IMPLICIT_INTERNAL_HPP
#define NEWTON_IMPLICIT_INTERNAL_HPP

#include "newton_implicit_solver.hpp"
#include "symbolic/fx/implicit_function_internal.hpp"
#include "symbolic/fx/nlp_solver.hpp"
#include "symbolic/fx/linear_solver.hpp"

namespace CasADi{

  /** \brief Internal class for NewtonImplicitInternal
   * 
      @copydoc ImplicitFunction_doc
   * */
class NewtonImplicitInternal : public ImplicitFunctionInternal {
  friend class NewtonImplicitSolver;
public:
  /** \brief  Constructor */
  explicit NewtonImplicitInternal();

  /** \brief  Clone */
  virtual NewtonImplicitInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit NewtonImplicitInternal(const FX& f);

  /** \brief  Destructor */
  virtual ~NewtonImplicitInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate(int nfdir, int nadir);

  /** \brief  Create a new ImplicitFunctionInternal */
  virtual ImplicitFunctionInternal* create(const FX& f) const { return new NewtonImplicitInternal(f);}
  
  protected:
    /// Maximum number of Newton iterations
    int max_iter_;
    
    /// Absolute tolerance that should be met on residual
    double abstol_;
    
    /// Absolute tolerance that should be met on step
    double abstolStep_;
};

} // namespace CasADi

#endif //NEWTON_IMPLICIT_INTERNAL_HPP

