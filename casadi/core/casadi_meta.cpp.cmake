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

#include <casadi/core/casadi_meta.hpp>

namespace casadi {
  const std::string CasadiMeta::version = "${PACKAGE_VERSION}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::git_revision = "${git_revision}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::git_describe = "${git_describe}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::feature_list = "${feature_list}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::build_type = "${CMAKE_BUILD_TYPE}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler_id = "${CMAKE_CXX_COMPILER_ID}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler = "${CMAKE_CXX_COMPILER}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler_flags ="${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::modules ="${CASADI_MODULES}";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::plugins ="${CASADI_PLUGINS}";  // NOLINT(whitespace/line_length)
} // namespace casadi
