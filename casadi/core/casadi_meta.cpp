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


#include "casadi_meta.hpp"
#include <casadi_meta.h>

namespace casadi {
  const char* CasadiMeta::version() { return casadi_version ;}

  const char* CasadiMeta::git_revision() { return casadi_git_revision ;}

  const char* CasadiMeta::git_describe() { return casadi_git_describe ;}

  const char* CasadiMeta::feature_list() { return casadi_feature_list ;}

  const char* CasadiMeta::build_type() { return casadi_build_type ;}

  const char* CasadiMeta::compiler_id() { return casadi_compiler_id ;}

  const char* CasadiMeta::compiler() { return casadi_compiler ;}

  const char* CasadiMeta::compiler_flags() { return casadi_compiler_flags ;}

  const char* CasadiMeta::modules() { return casadi_modules ;}

  const char* CasadiMeta::plugins() { return casadi_plugins ;}

  const char* CasadiMeta::install_prefix() { return casadi_install_prefix ;}

}  // namespace casadi
