/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef IP_INTERNAL_HPP
#define IP_INTERNAL_HPP

#include "ip_method.hpp"
#include "core/function/nlp_solver_internal.hpp"
#include "core/function/linear_solver.hpp"

namespace casadi{
    
class IPInternal : public NlpSolverInternal{

public:
  explicit IPInternal(const Function& F, const Function& G);
  virtual ~IPInternal();
  virtual IPInternal* clone() const{ return new IPInternal(*this);}
  
  virtual void init();
  virtual void evaluate(int nfdir, int nadir);
  
  /// Function which forms the KKT system
  enum KIn{K_x,K_t,K_NUM_IN};
  enum KOut{K_K,K_k,K_NUM_OUT};
  Function kfcn_;
  
  /// Linear solver for the KKT system
  LinearSolver linear_solver_;
};

} // namespace casadi

#endif //IP_INTERNAL_HPP
