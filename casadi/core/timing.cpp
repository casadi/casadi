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


#include "timing.hpp"

namespace casadi {

  timer getTimerTime() {
    timer ret;
    ret.user = clock();
    ret.real = getRealTime();
    return ret;
  }

  // ret = t1 - t0
  diffTime diffTimers(const timer t1, const timer t0) {
    diffTime ret;
    ret.user = (t1.user - t0.user)/CLOCKS_PER_SEC;
    ret.real = t1.real - t0.real;
    return ret;
  }

  // t += diff
  void timerPlusEq(diffTime & t, const diffTime diff) {
    t.user += diff.user;
    t.real += diff.real;
  }

  Dict diffToDict(const diffTime& diff) {
    Dict ret;
    // compatable names with the linux "time" utility
    ret["real"] = diff.real;
    ret["user"] = diff.user;
    return ret;
  }

} // namespace casadi
