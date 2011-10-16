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

#ifndef CASADI_OCP_HPP
#define CASADI_OCP_HPP

#include "casadi/casadi.hpp"

namespace CasADi{
  namespace OptimalControl{

/** \brief Symbolic representation of an optimal control problem (OCP) 
  Variables:
  t:     time
  x:     differential states
  xdot:  state derivatives
  z:     algebraic states
  y:     dependent variables
  u:     control signals
  p:     independent parameters
  
  Equations
  fully implicit DAE:       0 == dae(t,x,xdot,z,u,p)
  initial equations:        0 == ieq(t,x,xdot,z,u,p)
  explicit ODE:          xdot == ode(t,x,z,u,p)
  algebraic equations:      0 == alg(t,x,z,u,p)
  dependent equations:      y == dep(t,x,z,u,p)
  
  \author Joel Andersson 2011
*/
struct SymbolicOCP{
  
  /// Time
  SXMatrix t;
  
  /// Differential states
  SXMatrix x;

  /// State derivative
  SXMatrix xdot;

  /// Algebraic states
  SXMatrix z;
  
  /// Controls
  SXMatrix u;
      
  /// Free parameters
  SXMatrix p;
  
  /// Dependent variables
  SXMatrix y;

  /// Fully implicit DAE
  SXMatrix dae;

  /// Initial equations
  SXMatrix ieq;

  /// Explicit ODE
  SXMatrix ode;
  
  /// Algebraic equations
  SXMatrix alg;
  
  /// Dependent equations
  SXMatrix dep;
  
};

  } // namespace OptimalControl
} // namespace CasADi

#endif // CASADI_OCP_HPP


