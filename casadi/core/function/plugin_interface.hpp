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

#ifndef CASADI_PLUGIN_INTERFACE_HPP
#define CASADI_PLUGIN_INTERFACE_HPP

#include "../function/function_internal.hpp"
#include "../function/adaptor.hpp"
#include "../function/wrapper.hpp"

/// \cond INTERNAL

// For dynamic loading
#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#define NOMINMAX
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32

// Set default shared library prefix
#ifndef SHARED_LIBRARY_PREFIX
#define SHARED_LIBRARY_PREFIX "lib"
#endif // SHARED_LIBRARY_PREFIX

// Set default shared library suffix
#ifndef SHARED_LIBRARY_SUFFIX
#define SHARED_LIBRARY_SUFFIX ".so"
#endif // SHARED_LIBRARY_SUFFIX

#endif // WITH_DL

namespace casadi {

  typedef void (*AdaptorLoader)(const std::string& name);

  /** \brief Interface for accessing input and output data structures
      \author Joel Andersson
      \date 2013
  */
  template<class Derived>
  class CASADI_CORE_EXPORT PluginInterface {
    public:

    /// Fields
    struct Plugin{
      typename Derived::Creator creator;
      AdaptorLoader adaptorLoader;
      const char* name;
      const char* doc;
      int version;
    };

    // Plugin registration function
    typedef int (*RegFcn)(Plugin* plugin);

    /// Load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Register an integrator in the factory
    static void registerPlugin(RegFcn regfcn, const std::string& suffix = std::string());

    /// Load and get the creator function
    static Plugin& getPlugin(const std::string& name);

    /// Modifies the Derived object by setting Adaptor options when needed
    Derived* adaptor(const std::string &name_);

  };

  template<class Derived>
  void PluginInterface<Derived>::loadPlugin(const std::string& name_) {

    std::string::size_type dotpos = name_.find(".");

    std::string name;
    std::string suffix;

    if (dotpos == std::string::npos) {
      name = name_;
    } else {
      name = name_.substr(0, dotpos);
      suffix = name_.substr(dotpos+1);
    }

    // Issue warning aRegFcn reg;nd quick return if already loaded
    if (Derived::solvers_.find(name) != Derived::solvers_.end()) {
      casadi_warning("PluginInterface: Solver " + name + " is already in use. Ignored.");
    } else {

#ifndef WITH_DL
    casadi_error("WITH_DL option needed for dynamic loading");
#else // WITH_DL
    // Retrieve the registration function
    RegFcn reg;

    // Get the name of the shared library
    std::string lib = SHARED_LIBRARY_PREFIX "casadi_"
      + Derived::infix_ + "_" + name + SHARED_LIBRARY_SUFFIX;

    // Load the dll
    std::string regName = "casadi_register_" + Derived::infix_ + "_" + name;

    // Error string
    std::string errors = "PluginInterface::loadPlugin: Cannot load shared library:";
#ifdef _WIN32
    HINSTANCE handle = LoadLibrary(TEXT(lib.c_str()));
    if (!handle) {
      errors += "\n  Tried " + lib + ":\n    Error code (WIN32): " + STRING(GetLastError());

      #ifdef PLUGIN_EXTRA_SEARCH_PATH
      // Try the second search path
      std::string lib2 = PLUGIN_EXTRA_SEARCH_PATH "\\" + lib;
      handle = LoadLibrary(TEXT(lib2.c_str()));
      if (!handle) {
        errors += "\n  Tried: " + lib2 + ":\n    Error code (WIN32): " + STRING(GetLastError());
      }
      #endif // PLUGIN_EXTRA_SEARCH_PATH

      if (!handle) {
        // Try current directory
        std::string lib3 = ".\\" + lib;
        handle = LoadLibrary(TEXT(lib3.c_str()));
        if (!handle) {
          errors += "\n  Tried: " + lib3 + ":\n    Error code (WIN32): " + STRING(GetLastError());
        }
      }
    }
    casadi_assert_message(handle!=0, errors);

    reg = (RegFcn)GetProcAddress(handle, TEXT(regName.c_str()));
    casadi_assert_message(reg!=0, "PluginInterface::loadPlugin: no \"" + regName + "\" found");
#else // _WIN32
    void* handle = dlopen(lib.c_str(), RTLD_LAZY | RTLD_GLOBAL);
    if (!handle) {
      errors += "\n  Tried " + lib + ":\n    Error code: " + dlerror();

      #ifdef PLUGIN_EXTRA_SEARCH_PATH
      // Try the second search path
      lib = PLUGIN_EXTRA_SEARCH_PATH "/" + lib;
      handle = dlopen(lib.c_str(), RTLD_LAZY | RTLD_GLOBAL);
      if (!handle) {
        errors += "\n  Tried " + lib + ":\n    Error code: " + dlerror();
      }
      #endif // PLUGIN_EXTRA_SEARCH_PATH
    }
    casadi_assert_message(handle!=0, errors);

    // Reset error
    dlerror();

    // Load creator
    reg = (RegFcn)dlsym(handle, regName.c_str());
    casadi_assert_message(!dlerror(), "PluginInterface::loadPlugin: no \""+regName+"\" found");
#endif // _WIN32

    // Register the plugin
    registerPlugin(reg, suffix);
#endif // WITH_DL
    }

    if (suffix.size()>0) {
      Plugin &plugin = Derived::solvers_.find(name)->second;
      casadi_assert_message(plugin.adaptorLoader,
        "PluginInterface: could not find adaptor for " << name);
      plugin.adaptorLoader(suffix);
    }
  }

  template<class Derived>
  void PluginInterface<Derived>::registerPlugin(RegFcn regfcn, const std::string& suffix) {
    // Create a temporary struct
    Plugin plugin;
    plugin.adaptorLoader = 0;

    // Set the fields
    int flag = regfcn(&plugin);
    casadi_assert(flag==0);

    // Check if the solver name is in use
    typename std::map<std::string, Plugin>::iterator it=Derived::solvers_.find(plugin.name);
    casadi_assert_message(it==Derived::solvers_.end(),
                          "Solver " << plugin.name << " is already in use");

    // Add to list of solvers
    Derived::solvers_[plugin.name] = plugin;
  }

  template<class Derived>
  typename PluginInterface<Derived>::Plugin&
  PluginInterface<Derived>::getPlugin(const std::string& name_) {

    std::string::size_type dotpos = name_.find(".");

    std::string name;
    std::string suffix;

    if (dotpos == std::string::npos) {
      name = name_;
    } else {
      name = name_.substr(0, dotpos);
      suffix = name_.substr(dotpos+1);
    }

    // Check if the solver has been loaded
    typename std::map<std::string, Plugin>::iterator it=Derived::solvers_.find(name);

    // Load the solver if needed
    if (it==Derived::solvers_.end()) {
      loadPlugin(name);
      it=Derived::solvers_.find(name);
    }
    casadi_assert(it!=Derived::solvers_.end());
    return it->second;
  }

  template<class Derived>
  Derived* PluginInterface<Derived>::adaptor(const std::string &name_) {

    Derived* s = static_cast< Derived *>(this);
    AdaptorBase* dd = dynamic_cast< AdaptorBase *>(s);

    // Quick return if Derived instance is not an Adaptor
    if (!dd) {
      return s;
    }

    // Split name_ into: name + '.' + suffix
    std::string::size_type dotpos = name_.find(".");

    std::string name;
    std::string suffix;

    if (dotpos == std::string::npos) {
      name = name_;
    } else {
      name = name_.substr(0, dotpos);
      suffix = name_.substr(dotpos+1);
    }

    if (suffix.size()>0) {
      // If there is a suffix, use it to set the target Solver
      s->setOption(dd->prefix()+"_solver", suffix);
      s->setOption("target", false);
    } else {
      s->setOption("target", true);
    }

    return s;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_PLUGIN_INTERFACE_HPP
