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

#ifndef NLP_IMPLICIT_INTERNAL_HPP
#define NLP_IMPLICIT_INTERNAL_HPP

#include "symbolic/fx/implicit_function_internal.hpp"
#include "symbolic/fx/nlp_solver.hpp"

namespace CasADi{

  /** \brief Internal class for NLPImplicitInternal
   * 
      @copydoc ImplicitFunction_doc
   * */
class NLPImplicitInternal : public ImplicitFunctionInternal {
  friend class NLPImplicitSolver;
public:
  /** \brief  Constructor */
  explicit NLPImplicitInternal();

  /** \brief  Clone */
  virtual NLPImplicitInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit NLPImplicitInternal(const FX& f, int nrhs=1);

  /** \brief  Destructor */
  virtual ~NLPImplicitInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate(int nfdir, int nadir);

  protected:
    NLPSolver nlp_solver_;
};

} // namespace CasADi

#endif //NLP_IMPLICIT_INTERNAL_HPP

