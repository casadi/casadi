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


#include "ghc.hpp"
#include <ghc/filesystem.hpp>

namespace casadi {

    bool is_directory(const std::string& path) {
      return ghc::filesystem::is_directory(path);
    }

    bool create_directories(const std::string& path) {
      return ghc::filesystem::create_directories(path);
    }

    bool remove(const std::string& path) {
        return ghc::filesystem::remove(path);
    }

    casadi_int remove_all(const std::string& path) {
        return ghc::filesystem::remove_all(path);
    }

    std::string filename(const std::string& path) {
        return ghc::filesystem::path(path).filename().string();
    }

    std::string parent_path(const std::string& path) {
        return ghc::filesystem::path(path).parent_path().string();
    }

    bool has_parent_path(const std::string& path) {
        return ghc::filesystem::path(path).has_parent_path();
    }

    extern "C"
    int CASADI_FILESYSTEM_GHC_EXPORT
    casadi_register_filesystem_ghc(Filesystem::Plugin* plugin) {
        plugin->name = "ghc";
        plugin->doc = Ghc::meta_doc.c_str();
        plugin->version = CASADI_VERSION;
        plugin->exposed.is_directory = &is_directory;
        plugin->exposed.create_directories = &create_directories;
        plugin->exposed.remove = &remove;
        plugin->exposed.remove_all = &remove_all;
        plugin->exposed.filename = &filename;
        plugin->exposed.has_parent_path = &has_parent_path;
        plugin->exposed.parent_path = &parent_path;
        return 0;
    }

    extern "C"
    void CASADI_FILESYSTEM_GHC_EXPORT casadi_load_filesystem_ghc() {
        Filesystem::registerPlugin(casadi_register_filesystem_ghc);
    }


} // namespace casadi
