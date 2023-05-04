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


#include "timing.hpp"

namespace casadi {

  using namespace std::chrono;
  FStats::FStats() {
  }

  void FStats::reset() {
    n_call = 0;
    t_wall = 0;
    t_proc = 0;
  }

  void FStats::tic() {
    start_proc = std::clock();
    start_wall= high_resolution_clock::now();
  }

  void FStats::toc() {
    // First get the time points
    stop_proc = std::clock();
    stop_wall = high_resolution_clock::now();

    // Process them
    double proc = static_cast<double>(stop_proc - start_proc) / static_cast<double>(CLOCKS_PER_SEC);
    t_proc += proc;
    double wall = duration<double>(stop_wall - start_wall).count();
    t_wall += wall;
    n_call +=1;

  }

  void FStats::join(FStats& rhs) {
    t_proc += rhs.t_proc;
    t_wall += rhs.t_wall;
    n_call += rhs.n_call;
  }

  ScopedTiming::ScopedTiming(FStats& f) : f_(f) {
    f_.tic();
  }

  ScopedTiming::~ScopedTiming() {
    f_.toc();
  }

} // namespace casadi
