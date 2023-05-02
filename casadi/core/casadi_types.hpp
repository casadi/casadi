/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_TYPES_HPP
#define CASADI_TYPES_HPP

#ifndef CASADI_INT_TYPE
#define CASADI_INT_TYPE long long int
#endif // CASADI_INT_TYPE


#ifndef CASADI_C_H
// Avoid interaction problems with <casadi/casadi_c.h>
#ifndef CASADI_MEM_H
// Avoid interaction problems with <casadi/mem.h>
typedef CASADI_INT_TYPE casadi_int;
#endif
#endif

#ifndef SWIG  // Not in public API
typedef unsigned CASADI_INT_TYPE casadi_uint;
#endif // SWIG

// Convert to string
#define CASADI_STR1(x) #x
#define CASADI_STR(x) CASADI_STR1(x)
#define CASADI_INT_TYPE_STR CASADI_STR(CASADI_INT_TYPE)

typedef casadi_int casadi_index;

namespace casadi {
    /** Unified return status for solvers */
    enum UnifiedReturnStatus {
        SOLVER_RET_SUCCESS,
        SOLVER_RET_UNKNOWN,
        SOLVER_RET_LIMITED, // Out of time
        SOLVER_RET_NAN,
        SOLVER_RET_INFEASIBLE
    };
}


#endif // CASADI_TYPES_HPP
