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

#include "casadi_meta.hpp"

namespace casadi {
  const std::string CasadiMeta::version = "1.9.0+";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::git_revision = "6ae700a605fd7655e73b6a30c1d5c3b01bfc3642";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::git_describe = "1.9.0+258.6ae700a";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::feature_list = "\n * dynamic-loading , Compile with support for dynamic loading of generated functions (needed for ExternalFunction)\n * using-c++11 , Using C++11 features (improves efficiency and is required for some examples).\n * sundials-interface , Interface to the ODE/DAE integrator suite SUNDIALS.\n * csparse-interface , Interface to the sparse direct linear solver CSparse.\n * lapack-interface , Interface to LAPACK.\n * ipopt-interface , Interface to the NLP solver Ipopt.\n * snopt-interface , Interface to SNOPT.\n * qpoases-interface , Interface to the active-set QP solver qpOASES.\n * ooqp-interface , Interface to the QP solver OOQP (requires BLAS and HSL libraries).\n * sqic-interface , Interface to the QP solver SQIC.\n * slicot-interface , Interface to the controls library SLICOT.\n * dsdp-interface , Interface to the interior point SDP solver DSDP (requires BLAS and LAPACK).\n";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::build_type = "Release";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler_id = "GNU";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler = "/usr/bin/g++-4.7";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler_flags =" -std=gnu++11 -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-delete-non-virtual-dtor -fPIC -O3 -DNDEBUG";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::modules ="casadi_core;casadi_lapack_interface;casadi_sundials_interface;casadi_ipopt_interface;casadi_qpoases_interface;casadi_csparse_interface;casadi_dsdp_interface;casadi_ooqp_interface;casadi_sqic_interface;casadi_snopt_interface;casadi_slicot_interface;casadi_nonlinear_programming;casadi_convex_programming;casadi_integration;casadi_optimal_control;casadi_control";  // NOLINT(whitespace/line_length)
} // namespace casadi
