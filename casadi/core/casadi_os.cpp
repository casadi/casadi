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

#include "casadi_os.hpp"
#include "exception.hpp"
#include "global_options.hpp"

namespace casadi {

// http://stackoverflow.com/questions/303562/c-format-macro-inline-ostringstream
#define STRING(ITEMS) \
  ((dynamic_cast<std::ostringstream &>(std::ostringstream() \
    . seekp(0, std::ios_base::cur) << ITEMS)) . str())

char pathsep() {
    #ifdef _WIN32
        return ';';
    #else
        return ':';
    #endif
}
std::string filesep() {
    #ifdef _WIN32
        return "\\";
    #else
        return "/";
    #endif
}

std::vector<std::string> get_search_paths() {

    // Build up search paths;
    std::vector<std::string> search_paths;

    // Search path: global casadipath option
    std::stringstream casadipaths(GlobalOptions::getCasadiPath());
    std::string casadipath;
    while (std::getline(casadipaths, casadipath, pathsep())) {
        search_paths.push_back(casadipath);
    }

    // Search path: CASADIPATH env variable
    char* pLIBDIR;
    pLIBDIR = getenv("CASADIPATH");

    if (pLIBDIR!=nullptr) {
        std::stringstream casadipaths(pLIBDIR);
        std::string casadipath;
        while (std::getline(casadipaths, casadipath, pathsep())) {
        search_paths.push_back(casadipath);
        }
    }

    // Search path: bare
    search_paths.push_back("");

    // Search path : PLUGIN_EXTRA_SEARCH_PATH
    #ifdef PLUGIN_EXTRA_SEARCH_PATH
    search_paths.push_back(
        std::string("") + PLUGIN_EXTRA_SEARCH_PATH);
    #endif // PLUGIN_EXTRA_SEARCH_PATH

    // Search path : current directory
    search_paths.push_back(".");

    return search_paths;
}

#ifdef WITH_DL

handle_t open_shared_library(const std::string& lib, const std::vector<std::string> &search_paths,
    const std::string& caller, bool global) {
        std::string resultpath;
        return open_shared_library(lib, search_paths, resultpath, caller, global);
}

int close_shared_library(handle_t handle) {
    #ifdef _WIN32
        return !FreeLibrary(handle);
    #else // _WIN32
        return dlclose(handle);
    #endif // _WIN32
}

handle_t open_shared_library(const std::string& lib, const std::vector<std::string> &search_paths,
        std::string& resultpath, const std::string& caller, bool global) {
    // Alocate a handle
    handle_t handle = 0;

    // Alocate a handle pointer
    #ifndef _WIN32
        int flag;
        if (global) {
            flag = RTLD_NOW | RTLD_GLOBAL;
        } else {
            flag = RTLD_LAZY | RTLD_LOCAL;
        }
    #ifdef WITH_DEEPBIND
    #ifndef __APPLE__
        flag |= RTLD_DEEPBIND;
    #endif
    #endif
    #endif


    // Prepare error string
    std::stringstream errors;
    errors << caller << ": Cannot load shared library '"
           << lib << "': " << std::endl;
    errors << "   (\n"
           << "    Searched directories: 1. casadipath from GlobalOptions\n"
           << "                          2. CASADIPATH env var\n"
           << "                          3. PATH env var (Windows)\n"
           << "                          4. LD_LIBRARY_PATH env var (Linux)\n"
           << "                          5. DYLD_LIBRARY_PATH env var (osx)\n"
           << "    A library may be 'not found' even if the file exists:\n"
           << "          * library is not compatible (different compiler/bitness)\n"
           << "          * the dependencies are not found\n"
           << "   )";

    std::string searchpath;

    // Try getting a handle
    for (casadi_int i=0;i<search_paths.size();++i) {
      searchpath = search_paths[i];
#ifdef _WIN32
      SetDllDirectory(TEXT(searchpath.c_str()));
      handle = LoadLibrary(TEXT(lib.c_str()));
      SetDllDirectory(NULL);
#else // _WIN32
      std::string libname = searchpath.empty() ? lib : searchpath + filesep() + lib;
      handle = dlopen(libname.c_str(), flag);
#endif // _WIN32
      if (handle) {
        resultpath = searchpath;
        break;
      } else {
        errors << std::endl << "  Tried '" << searchpath << "' :";
#ifdef _WIN32
        errors << std::endl << "    Error code (WIN32): " << STRING(GetLastError());
#else // _WIN32
        errors << std::endl << "    Error code: " << dlerror();
#endif // _WIN32
      }
    }

    casadi_assert(handle!=nullptr, errors.str());

    return handle;
}

#endif // WITH_DL

} // namespace casadi
