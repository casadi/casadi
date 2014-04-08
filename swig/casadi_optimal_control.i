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
 
%module casadi_optimal_control

%include "common.i"

%import "casadi_symbolic.i"

#define CASADI_OPTIMAL_CONTROL_EXPORT

%{
#include "casadi/optimal_control/variable.hpp"
#include "casadi/optimal_control/symbolic_ocp.hpp"
#include "casadi/optimal_control/direct_single_shooting.hpp"
#include "casadi/optimal_control/direct_multiple_shooting.hpp"
#include "casadi/optimal_control/direct_collocation.hpp"
%}

%include "casadi/optimal_control/variable.hpp"
%include "casadi/optimal_control/symbolic_ocp.hpp"
%include "casadi/optimal_control/direct_single_shooting.hpp"
%include "casadi/optimal_control/direct_multiple_shooting.hpp"
%include "casadi/optimal_control/direct_collocation.hpp"

#ifdef SWIGPYTHON
%pythoncode %{
  class VariableStruct(object):
    """Structure for browsing through a variable tree."""
    def __repr__(self):
      return repr(self.__dict__)
%}

#endif // SWIGPYTHON

VECTOR_REPR(casadi::Variable)


