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


#include "casadi_options.hpp"
#include "casadi_exception.hpp"

namespace casadi {

  bool CasadiOptions::catch_errors_swig = true;
  bool CasadiOptions::simplification_on_the_fly = true;
  bool CasadiOptions::profiling = false;
  std::ofstream CasadiOptions::profilingLog;
  bool CasadiOptions::profilingBinary = true;
  bool CasadiOptions::purgeSeeds = false;
  bool CasadiOptions::allowed_internal_api = false;

  void CasadiOptions::startProfiling(const std::string &filename) {
    profilingLog.open(filename.c_str(), std::ofstream::out);
    if (profilingLog.is_open()) {
      profiling = true;
    } else {
      casadi_error("Did not manage to open file " << filename << " for logging.");
    }
  }

  void CasadiOptions::stopProfiling() {
    if (profiling) {
      profilingLog.close();
    }
    profiling = false;
  }

} // namespace casadi
