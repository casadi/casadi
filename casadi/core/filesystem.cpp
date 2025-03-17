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


#include "filesystem.hpp"

namespace casadi {

std::map<std::string, Filesystem::Plugin> Filesystem::solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
std::mutex Filesystem::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

const std::string Filesystem::infix_ = "filesystem";

bool Filesystem::is_directory(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.is_directory(path);
}

bool Filesystem::remove(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.remove(path);
}

casadi_int Filesystem::remove_all(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.remove_all(path);
}

std::string Filesystem::filename(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.filename(path);
}

bool Filesystem::is_enabled() {
  return Filesystem::has_plugin("ghc");
}

void Filesystem::assert_enabled() {
  casadi_assert(Filesystem::is_enabled(),
  "This action requires advanced filesystem access. Compile CasADi with WITH_GC=ON.");
}

} // namespace casadi
