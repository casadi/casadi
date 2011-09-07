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

#ifndef LIFTOPT_INTERNAL_HPP
#define LIFTOPT_INTERNAL_HPP

#include <casadi/fx/mx_function.hpp>
#include <casadi/fx/nlp_solver_internal.hpp>
#include "liftopt_solver.hpp"

#include <liftopt.hpp>
#include <casadi/fx/mx_function.hpp>
#include <casadi/fx/sx_function.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/stl_vector_tools.hpp>
#include "sonic++.h"

namespace CasADi{
  namespace Interfaces{

/**

@copydoc NLPSolver_doc
*/
class LiftoptInternal : public NLPSolverInternal{
  public:
    LiftoptInternal(const MXFunction& fcn);
    virtual ~LiftoptInternal();
    virtual LiftoptInternal* clone() const{ return new LiftoptInternal(*this);}
    virtual void init();
    virtual void evaluate(int nfdir, int nadir);

    std::vector<double> nodeInit;
    
    MXFunction fcn_;
    liftopt::ILifter *problem_;
    liftopt::IOptimizer *opt_ ;
    liftopt::IQPSolver* qp_solver_;
    liftopt::DVec uInit_;
    liftopt::DVec lambdaInit_;
    liftopt::DVec loCtrlBounds_;
    liftopt::DVec upCtrlBounds_;
    liftopt::DVec nodeInit_;
    int m_nObjRes_, m_nCtrls_, m_nEq_, m_nIneq_;
};

  } // namespace Interfaces
} // namespace CasADi

#endif // LIFTOPT_INTERNAL_HPP
