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

#ifndef SIMPLE_HOMOTOPY_NLP_INTERNAL_HPP
#define SIMPLE_HOMOTOPY_NLP_INTERNAL_HPP

#include "simple_homotopy_nlp_solver.hpp"
#include "symbolic/function/homotopy_nlp_internal.hpp"
#include "symbolic/function/stabilized_qp_solver.hpp"

/// \cond INTERNAL
namespace CasADi{
    
class SimpleHomotopyNLPInternal : public HomotopyNLPInternal{

public:
  explicit SimpleHomotopyNLPInternal(const Function& hnlp);
  virtual ~SimpleHomotopyNLPInternal();
  virtual SimpleHomotopyNLPInternal* clone() const{ return new SimpleHomotopyNLPInternal(*this);}
  
  virtual void init();
  virtual void evaluate();
  
  NLPSolver nlpsolver_;
  
  /// Take this many steps to go from tau=0 to tau=1
  int num_steps_;
  
};
/// \endcond
} // namespace CasADi

#endif //SIMPLE_HOMOTOPY_NLP_INTERNAL_HPP
