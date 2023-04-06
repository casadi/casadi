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


#ifndef CASADI_TIMING_HPP
#define CASADI_TIMING_HPP

#include "generic_type.hpp"

#include <chrono>
#include <ctime>

namespace casadi {
  /// \cond INTERNAL

  /**
  Timer class


  FStats hack;
  hack.tic();
  ....
  hack.toc();

  */
  class CASADI_EXPORT FStats {
    private:
      /// Time point used for wall time computation
      std::chrono::time_point<std::chrono::high_resolution_clock> start_wall;

      /// Time point used for proc time computation
      std::clock_t start_proc;

      /// Time point used for wall time computation
      std::chrono::time_point<std::chrono::high_resolution_clock> stop_wall;

      /// Time point used for proc time computation
      std::clock_t stop_proc;

    public:
      /// Constructor
      FStats();

      /// Reset the statistics
      void reset();

      /// Start timing
      void tic();

      /// Stop timing
      void toc();

      /// Accumulated number of calls since last reset
      casadi_int n_call = 0;

      /// Accumulated wall time [s] since last reset
      double t_wall = 0;

      /// Accumulated proc time [s] since last reset
      double t_proc = 0;

      void join(FStats& rhs);

  };

  class CASADI_EXPORT ScopedTiming {
    public:
      ScopedTiming(FStats& f);
      ~ScopedTiming();
    private:
      FStats& f_;
  };

/// \endcond
} // namespace casadi

#endif // CASADI_TIMING_HPP
