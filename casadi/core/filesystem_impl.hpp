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


#ifndef CASADI_FILESYSTEM_IMPL_HPP
#define CASADI_FILESYSTEM_IMPL_HPP

#include "filesystem.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL
namespace casadi {

  class CASADI_EXPORT
  Filesystem : public PluginInterface<Filesystem> {
  public:
    typedef bool (* IsDirectory)(const std::string& path);
    typedef bool (* CreateDirectories)(const std::string& path);
    typedef bool (* Remove)(const std::string& path);
    typedef casadi_int (* RemoveAll)(const std::string& path);
    typedef std::string (* Filename)(const std::string& path);
    typedef  bool (* HasParentPath)(const std::string& path);
    typedef std::string (* ParentPath)(const std::string& path);

    // Creator function for internal class
    typedef Filesystem* (*Creator)();

    static const std::string meta_doc;

    // No static functions exposed
    struct Exposed{
      IsDirectory is_directory;
      CreateDirectories create_directories;
      Remove remove;
      RemoveAll remove_all;
      Filename filename;
      HasParentPath has_parent_path;
      ParentPath parent_path;
    };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    static void assert_enabled();
    static bool is_directory(const std::string& path);
    static bool remove(const std::string& path);
    static casadi_int remove_all(const std::string& path);
    static std::string filename(const std::string& path);
    static bool is_enabled();
    static bool has_parent_path(const std::string& path);
    static std::string parent_path(const std::string& path);
    static bool ensure_directory_exists(const std::string& path);
    static bool create_directories(const std::string& path);

    static void open(std::ofstream&, const std::string& path,
        std::ios_base::openmode mode = std::ios_base::out);

    static std::ofstream* ofstream_ptr(const std::string& path,
        std::ios_base::openmode mode = std::ios_base::out);

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /// Infix
    static const std::string infix_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_FILESYSTEM_IMPL_HPP
