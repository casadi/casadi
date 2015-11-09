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
#include <ctime>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

namespace casadi {
  Timer getTimerTime() {
    Timer ret;
    ret.user = clock();
#ifdef _WIN32
    FILETIME tm;
    ULONGLONG t;
    GetSystemTimePreciseAsFileTime(&tm);
    t = (static_cast<ULONGLONG>(tm.dwHighDateTime) << 32) | (ULONGLONG)tm.dwLowDateTime;
    ret.real = static_cast<double>(t) / 10000000.0;
#else
    struct timeval tm;
    gettimeofday(&tm, NULL);
    ret.real = tm.tv_sec + tm.tv_usec/1000000.0;
#endif
    return ret;
  }

  // ret = t1 - t0
  DiffTime diffTimers(const Timer t1, const Timer t0) {
    DiffTime ret;
    ret.user = (t1.user - t0.user)/CLOCKS_PER_SEC;
    ret.real = t1.real - t0.real;
    return ret;
  }

  // t += diff
  void timerPlusEq(DiffTime & t, const DiffTime diff) {
    t.user += diff.user;
    t.real += diff.real;
  }

  Dict diffToDict(const DiffTime& diff) {
    Dict ret;
    // compatable names with the linux "time" utility
    ret["real"] = diff.real;
    ret["user"] = diff.user;
    return ret;
  }
} // namespace casadi
