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


#include "shell_compiler.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/casadi_meta.hpp"
#include <fstream>
#include <stdlib.h>
#include <dlfcn.h>
#include <cstdlib>
#include <unistd.h>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_COMPILER_SHELL_EXPORT
  casadi_register_compiler_shell(CompilerInternal::Plugin* plugin) {
    plugin->creator = ShellCompiler::creator;
    plugin->name = "shell";
    plugin->doc = ShellCompiler::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_COMPILER_SHELL_EXPORT casadi_load_compiler_shell() {
    CompilerInternal::registerPlugin(casadi_register_compiler_shell);
  }

  ShellCompiler* ShellCompiler::clone() const {
    // Return a deep copy
    ShellCompiler* node = new ShellCompiler(name_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  ShellCompiler::ShellCompiler(const std::string& name) :
    CompilerInternal(name) {
    addOption("compiler", OT_STRING, "gcc", "Compiler command");
    addOption("compiler_setup", OT_STRING, "-fPIC -shared", "Compiler setup command");
    addOption("flags", OT_STRINGVECTOR, GenericType(),
      "Compile flags for the JIT compiler. Default: None");
  }

  ShellCompiler::~ShellCompiler() {
    // Unload
    if (handle_) dlclose(handle_);

    // Delete the temporary file
    std::string rmcmd = "rm " + bin_name_;
    if (system(rmcmd.c_str())) {
      casadi_warning("Failed to delete temporary file:" + bin_name_);
    }
  }

  void ShellCompiler::init() {
    // Initialize the base classes
    CompilerInternal::init();

    // Read options
    string compiler = getOption("compiler").toString();
    string compiler_setup = getOption("compiler_setup").toString();
    vector<string> flags;
    if (hasSetOption("flags")) flags = getOption("flags");

    // Construct the compiler command
    stringstream cmd;
    cmd << compiler << " " << compiler_setup;
    for (vector<string>::const_iterator i=flags.begin(); i!=flags.end(); ++i) {
      cmd << " " << *i;
    }

    // C/C++ source file
    cmd << " " << name_;

    // Temporary file
#ifdef HAVE_MKSTEMPS
    // Preferred solution
    char bin_name[] = "tmp_casadi_compiler_shell_XXXXXX.so";
    if (mkstemps(bin_name, 3) == -1) {
      casadi_error("Failed to create a temporary file name");
    }
    bin_name_ = bin_name;
#else
    // Fallback, may result in deprecation warnings
    char* bin_name = tempnam(0, "tmp_casadi_compiler_shell_");
    bin_name_ = bin_name;
    free(bin_name);
#endif

    cmd << " -o " << bin_name_;

    // Compile into a shared library
    if (system(cmd.str().c_str())) {
      casadi_error("Compilation failed. Tried \"" + cmd.str() + "\"");
    }

    // Load shared library
    handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open function: "
                          << bin_name_ << ". error code: "<< dlerror());
    // reset error
    dlerror();
  }

  void* ShellCompiler::getFunction(const std::string& symname) {
    void* ret;
    ret = reinterpret_cast<void*>(dlsym(handle_, symname.c_str()));
    if (dlerror()) {
      ret=0;
      dlerror(); // Reset error flags
    }
    return ret;
  }

} // namespace casadi
