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

#include "function_internal.hpp"
#include "global_options.hpp"
#include "serializing_stream.hpp"
#include "casadi_os.hpp"
#include <casadi/config.h>

#include <stdlib.h>

/// \cond INTERNAL

namespace casadi {
  // Avoid segmentation faults when exposed function not implemented
  template<typename T>
  T check_exposed(T t) {
    casadi_assert(t!=0, "Static function not implemented for plugin");
    return t;
  }

  typedef ProtoFunction* (*Deserialize)(DeserializingStream&);

  /** \brief Interface for accessing input and output data structures

      \author Joel Andersson
      \date 2013

      \identifier{rp} */
  template<class Derived>
  class PluginInterface {
    public:


    /// Fields
    struct Plugin{
      typename Derived::Creator creator;
      const char* name;
      const char* doc;
      int version;
      typename Derived::Exposed exposed;
      const Options* options;
      Deserialize deserialize;
      // Constructor
      Plugin() : creator(nullptr), name(nullptr), doc(nullptr), version(0),
                  options(nullptr), deserialize(nullptr) {}
    };

    // Plugin registration function
    typedef int (*RegFcn)(Plugin* plugin);

    /// Check if a plugin is available or can be loaded
    static bool has_plugin(const std::string& pname, bool verbose=false);

    /// Get the plugin options
    static const Options& plugin_options(const std::string& pname);

    /// Get the plugin deserialize_map
    static Deserialize plugin_deserialize(const std::string& pname);

    /// Instantiate a Plugin struct from a factory function
    static Plugin pluginFromRegFcn(RegFcn regfcn);

    /// Load a plugin dynamically
    static Plugin load_plugin(const std::string& pname, bool register_plugin=true);

    /// Load a library dynamically
    static handle_t load_library(const std::string& libname, std::string& resultpath,
      bool global);

    /// Register an integrator in the factory
    static void registerPlugin(const Plugin& plugin);

    /// Register an integrator in the factory
    static void registerPlugin(RegFcn regfcn);

    /// Load and get the creator function
    static Plugin& getPlugin(const std::string& pname);

    // Create solver instance
    template<class Problem>
      static Derived* instantiate(const std::string& fname,
                                        const std::string& pname, Problem problem);
    // Get name of the plugin
    virtual const char* plugin_name() const = 0;

    /** \brief Serialize type information

        \identifier{rq} */
    void serialize_type(SerializingStream& s) const {
      s.pack("PluginInterface::plugin_name", std::string(plugin_name()));
    }

    /** \brief Deserialize with type disambiguation

        \identifier{rr} */
    static ProtoFunction* deserialize(DeserializingStream& s) {
      std::string class_name, plugin_name;
      s.unpack("PluginInterface::plugin_name", plugin_name);
      Deserialize deserialize = plugin_deserialize(plugin_name);
      return deserialize(s);
    }

  };

  template<class Derived>
  bool PluginInterface<Derived>::has_plugin(const std::string& pname, bool verbose) {

    // Quick return if available
    if (Derived::solvers_.find(pname) != Derived::solvers_.end()) {
      return true;
    }

    // Try loading the plugin
    try {
      (void)load_plugin(pname, false);
      return true;
    } catch (CasadiException& ex) {
      if (verbose) {
        casadi_warning(ex.what());
      }
      return false;
    }
  }

  template<class Derived>
  const Options& PluginInterface<Derived>::plugin_options(const std::string& pname) {
    const Options *op = getPlugin(pname).options;
    casadi_assert(op!=nullptr, "Plugin \"" + pname + "\" does not support options");
    return *op;
  }

  template<class Derived>
  Deserialize PluginInterface<Derived>::plugin_deserialize(const std::string& pname) {
    Deserialize m = getPlugin(pname).deserialize;
    casadi_assert(m, "Plugin \"" + pname + "\" does not support deserialize");
    return m;
  }

  template<class Derived>
  typename PluginInterface<Derived>::Plugin
      PluginInterface<Derived>::pluginFromRegFcn(RegFcn regfcn) {
    // Create a temporary struct
    Plugin plugin;

    // Set the fields
    int flag = regfcn(&plugin);
    casadi_assert(flag==0, "Registration of plugin failed.");

    return plugin;
  }


  template<class Derived>
  handle_t PluginInterface<Derived>::load_library(const std::string& libname,
    std::string& resultpath, bool global) {

#ifndef WITH_DL
    casadi_error("WITH_DL option needed for dynamic loading");
#else // WITH_DL

    // Get the name of the shared library
    std::string lib = SHARED_LIBRARY_PREFIX + libname + SHARED_LIBRARY_SUFFIX;

    // Build up search paths;
    std::vector<std::string> search_paths = get_search_paths();
    return open_shared_library(lib, search_paths, resultpath,
      "PluginInterface::load_plugin", global);

#endif // WITH_DL
  }

  template<class Derived>
  typename PluginInterface<Derived>::Plugin
      PluginInterface<Derived>::load_plugin(const std::string& pname, bool register_plugin) {
    // Issue warning and quick return if already loaded
    if (Derived::solvers_.find(pname) != Derived::solvers_.end()) {
      casadi_warning("PluginInterface: Solver " + pname + " is already in use. Ignored.");
      return Plugin();
    }


#ifndef WITH_DL
    casadi_error("WITH_DL option needed for dynamic loading");
#else // WITH_DL
    // Retrieve the registration function
    RegFcn reg;

    // Load the dll
    std::string regName = "casadi_register_" + Derived::infix_ + "_" + pname;

    std::string searchpath;
    handle_t handle = load_library("casadi_" + Derived::infix_ + "_" + pname, searchpath,
      false);

#ifdef _WIN32

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
    reg = reinterpret_cast<RegFcn>(GetProcAddress(handle, TEXT(regName.c_str())));
#if __GNUC__
#pragma GCC diagnostic pop
#endif

#else // _WIN32
    // Reset error
    dlerror();

    // Load creator
    reg = reinterpret_cast<RegFcn>(dlsym(handle, regName.c_str()));
#endif // _WIN32
    casadi_assert(reg!=nullptr,
      "PluginInterface::load_plugin: no \"" + regName + "\" found in " + searchpath + ".");

    // Create a temporary struct
    Plugin plugin = pluginFromRegFcn(reg);
    // Register the plugin
    if (register_plugin) {
      registerPlugin(plugin);
    }

    return plugin;

#endif // WITH_DL
  }

  template<class Derived>
  void PluginInterface<Derived>::registerPlugin(RegFcn regfcn) {
    registerPlugin(pluginFromRegFcn(regfcn));
  }

  template<class Derived>
  void PluginInterface<Derived>::registerPlugin(const Plugin& plugin) {

    // Check if the solver name is in use
    typename std::map<std::string, Plugin>::iterator it=Derived::solvers_.find(plugin.name);
    casadi_assert(it==Derived::solvers_.end(),
      "Solver " + str(plugin.name) + " is already in use");

    // Add to list of solvers
    Derived::solvers_[plugin.name] = plugin;
  }

  template<class Derived>
  typename PluginInterface<Derived>::Plugin&
  PluginInterface<Derived>::getPlugin(const std::string& pname) {

    // Check if the solver has been loaded
    auto it=Derived::solvers_.find(pname);

    // Load the solver if needed
    if (it==Derived::solvers_.end()) {
      load_plugin(pname);
      it=Derived::solvers_.find(pname);
    }
    casadi_assert_dev(it!=Derived::solvers_.end());
    return it->second;
  }

  template<class Derived>
  template<class Problem>
  Derived* PluginInterface<Derived>::
  instantiate(const std::string& fname,
              const std::string& pname, Problem problem) {

    // Assert the plugin exists (needed for adaptors)
    if (!has_plugin(pname, true)) {
      casadi_error("Plugin '" + pname + "' is not found.");
    }
    return getPlugin(pname).creator(fname, problem);
  }

} // namespace casadi

/// \endcond

#endif // CASADI_PLUGIN_INTERFACE_HPP
