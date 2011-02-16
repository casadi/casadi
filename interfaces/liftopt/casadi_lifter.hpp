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

#ifndef CASADI_LIFTER_HPP
#define CASADI_LIFTER_HPP

#include "liftopt_internal.hpp"

namespace CasADi{
  namespace Interfaces{

class CasadiLifter: public liftopt::ILifter{
  public:
    CasadiLifter(LiftoptInternal* interface): liftopt::ILifter( "casadi_interface" ), interface_(interface) {
      m_nObjRes = interface_->m_nObjRes_;
      m_nCtrls  = interface_->m_nCtrls_;
      m_nEq     = interface_->m_nEq_;
      m_nIneq   = interface_->m_nIneq_;
      m_presentProblemFormulations.set ( Problem_Residual );
      m_presentProblemFormulations.set ( Problem_LagGrad );
    };

    virtual ~CasadiLifter(){};
    virtual long evalUserFcn ( liftopt::TLifterArgs<double>& args );
    static void liftfun(double *v, int n, void *user_data);
    LiftoptInternal* interface_;
};

  } // namespace Interfaces
} // namespace CasADi


#endif // CASADI_LIFTER_HPP

