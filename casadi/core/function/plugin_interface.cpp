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

#include "plugin_interface.hpp"

namespace casadi {

  std::vector<std::string> getPluginSearchPaths() {
    // Build up search paths;
    std::vector<std::string> search_paths;

    #ifdef _WIN32
    char pathsep = ';';
    const std::string filesep("\\");
    #else
    char pathsep = ':';
    const std::string filesep("/");
    #endif

    // Search path: global casadipath option
    std::stringstream casadipaths(CasadiOptions::getCasadiPath());
    std::string casadipath;
    while (std::getline(casadipaths, casadipath, pathsep)) {
      search_paths.push_back(casadipath);
    }

    // Search path: CASADIPATH env variable
    char* pLIBDIR;
    pLIBDIR = getenv("CASADIPATH");

    if (pLIBDIR!=0) {
      std::stringstream casadipaths(pLIBDIR);
      std::string casadipath;
      while (std::getline(casadipaths, casadipath, pathsep)) {
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

} // namespace casadi
