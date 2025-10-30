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


#include "global_options.hpp"
#include "exception.hpp"
#include "filesystem_impl.hpp"
#include "casadi_meta.hpp"

namespace casadi {

  bool GlobalOptions::simplification_on_the_fly = true;
  bool GlobalOptions::hierarchical_sparsity = true;

  std::string GlobalOptions::casadipath;
  std::string GlobalOptions::casadi_include_path;
  std::string GlobalOptions::casadi_executable_path = CasadiMeta::bin_prefix();

  casadi_int GlobalOptions::max_num_dir = 64;

  // By default, use zero-based indexing
  casadi_int GlobalOptions::start_index = 0;

  std::string GlobalOptions::temp_work_dir = "./";

  bool GlobalOptions::julia_initialized = false;

  casadi_int GlobalOptions::copy_elision_min_size = 8;

  void GlobalOptions::setTempWorkDir(const std::string& dir) {
    casadi_assert(!dir.empty(), "Temporary working directory must be non-empty.");
    temp_work_dir = Filesystem::ensure_trailing_slash(dir);
  }
  casadi_int GlobalOptions::vector_width_real = 1;

  casadi_int GlobalOptions::byte_width_real = 8;


  bool GlobalOptions::feature_ve = true;
} // namespace casadi
