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


#ifndef CASADI_TIMING_HPP
#define CASADI_TIMING_HPP

#include "profiling.hpp"
#include "generic_type.hpp"

namespace casadi {
/// \cond INTERNAL

  typedef struct {
    double proc;
    double wall;
  } timer;

  typedef struct {
    double proc;
    double wall;
  } diffTime;

  CASADI_EXPORT timer getTimerTime(void);
  // ret = t1 - t0
  CASADI_EXPORT diffTime diffTimers(const timer t1, const timer t0);
  // t += diff
  CASADI_EXPORT void timerPlusEq(diffTime & t, const diffTime diff);

  CASADI_EXPORT Dict diffToDict(const diffTime& diff);

/// \endcond
} // namespace casadi

#endif // CASADI_TIMING_HPP
