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

#ifndef MULTIPLE_SHOOTING_INTERNAL_HPP
#define MULTIPLE_SHOOTING_INTERNAL_HPP

#include "multiple_shooting.hpp"
#include "../casadi/fx/ocp_solver_internal.hpp"

#include "../casadi/fx/parallelizer.hpp"
#include "../casadi/fx/c_function.hpp"
#include "../casadi/fx/mx_function.hpp"
#include "../casadi/fx/sx_function.hpp"

namespace CasADi{
  namespace OptimalControl{
  
class MultipleShootingInternal : public OCPSolverInternal{
  friend class MultipleShooting;
  
  public:
    // Constructor
    MultipleShootingInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn);
    
    // Initialize
    virtual void init();
    
    // NLP constraint
    static void constraint_wrapper(CFunction &f, int fsens_order, int asens_order, void* user_data);
    void constraint(CFunction &f, int fsens_order, int asens_order);

    // Jacobian of the NLP constraint
    static void jacobian_wrapper(CFunction &f, int fsens_order, int asens_order, void* user_data);
    void jacobian(CFunction &f, int fsens_order, int asens_order);

    // Solve the OCP
    virtual void evaluate(int fsens_order, int asens_order);
    
    FX ffcn2_;
  protected:
        
    // NLP objective function
    MXFunction F_;
    
    // NLP constraint function
    MXFunction G_;

    // Jacobian of the NLP constraints
    CFunction J_;

    // NLP solver
    NLPSolver nlp_solver_;

    // Parallel constraint function evaluation
    Parallelizer pF_;

    // Parallel evaluation of the Jacobian blocks
    Parallelizer pJX_,pJP_;
    
    // Mapping from the Jacobian blocks to the sparse Jacobian
    SXFunction J_mapping_;
};
                        
                        
  } // namespace OptimalControl
} // namespace CasADi


#endif // MULTIPLE_SHOOTING_INTERNAL_HPP
