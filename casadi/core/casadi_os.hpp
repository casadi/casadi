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

#ifndef CASADI_CASADI_OS_HPP
#define CASADI_CASADI_OS_HPP

#include <vector>
#include <string>

/// \cond INTERNAL


// For dynamic loading
#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#endif
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32

// Set default shared library prefix
#ifndef SHARED_LIBRARY_PREFIX
#define SHARED_LIBRARY_PREFIX CASADI_SHARED_LIBRARY_PREFIX
#endif // SHARED_LIBRARY_PREFIX


// Set default shared library suffix
#ifndef SHARED_LIBRARY_SUFFIX
#define SHARED_LIBRARY_SUFFIX CASADI_SHARED_LIBRARY_SUFFIX
#endif // SHARED_LIBRARY_SUFFIX
#endif // WITH_DL

namespace casadi {

std::vector<std::string> get_search_paths();

/* \brief Get the file separator (: or ;)
 */
char pathsep();
/* \brief Get the file separator (/ or \)
 */
std::string filesep();

// For dynamic loading
#ifdef WITH_DL

#ifdef _WIN32
    typedef HINSTANCE handle_t;
#else // _WIN32
    typedef void* handle_t;
#endif

handle_t open_shared_library(const std::string& lib, const std::vector<std::string> &search_paths,
    std::string &resultpath,
    const std::string& caller, bool global=false);

handle_t open_shared_library(const std::string& lib, const std::vector<std::string> &search_paths,
    const std::string& caller, bool global=false);

/** \brief Close shared library 
 * \return 0 if successful
 */
int close_shared_library(handle_t handle);

#endif // WITH_DL



} // namespace casadi

/// \endcond

#endif // CASADI_SHARED_LIBRARIES_HPP
