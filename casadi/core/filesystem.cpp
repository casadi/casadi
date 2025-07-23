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


#include "filesystem_impl.hpp"

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

bool Filesystem::has_parent_path(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.has_parent_path(path);
}

std::string Filesystem::parent_path(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.parent_path(path);
}

bool Filesystem::create_directories(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.create_directories(path);
}

std::vector<std::string> Filesystem::iterate_directory_names(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.iterate_directory_names(path);
}

std::string Filesystem::filename(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.filename(path);
}

std::string Filesystem::absolute(const std::string& path) {
  assert_enabled();
  return Filesystem::getPlugin("ghc").exposed.absolute(path);
}

bool Filesystem::is_enabled() {
  return Filesystem::has_plugin("ghc");
}

void Filesystem::assert_enabled() {
  casadi_assert(Filesystem::is_enabled(),
  "This action requires advanced filesystem access. Compile CasADi with WITH_GHC_FILESYSTEM=ON.");
}

bool has_filesystem(const std::string& name) {
  return Filesystem::has_plugin(name);
}

void load_filesystem(const std::string& name) {
  Filesystem::load_plugin(name);
}

std::string doc_filesystem(const std::string& name) {
  return Filesystem::getPlugin(name).doc;
}


bool Filesystem::ensure_directory_exists(const std::string& filename) {
  if (has_parent_path(filename)) {
    std::string dir = parent_path(filename);
    if (!is_directory(dir)) {
      return create_directories(dir);
    }
  }
  return true;
}

std::unique_ptr<std::ostream> Filesystem::ofstream_ptr(const std::string& path,
  std::ios_base::openmode mode) {
  if (is_enabled()) {
    casadi_assert(ensure_directory_exists(path),
     "Unable to create the required directory for '" + path + "'.");
  }
  auto stream_ptr = ofstream_compat(path, mode);
  if (is_enabled()) {
    casadi_assert(stream_ptr,
      "Error opening stream '" + path + "'.");
  } else {
    casadi_assert(stream_ptr,
      "Error opening stream '" + path + "'. "
      "Does the directory exits? "
      "Note that CasADi needs to be compiled with WITH_GHC_FILESYSTEM=ON "
      "for directories to be automatically created");
  }
  return stream_ptr;
}

std::unique_ptr<std::istream> Filesystem::ifstream_ptr(const std::string& path,
  std::ios_base::openmode mode, bool fail) {
  auto ret = ifstream_compat(path, mode);
  if (fail && !ret) {
    casadi_error("Could not open '" + path + "'.");
  }
  return ret;
}

bool Filesystem::exists(const std::string& path) {
  auto ret = ifstream_compat(path);
  return static_cast<bool>(ret);
}

} // namespace casadi
