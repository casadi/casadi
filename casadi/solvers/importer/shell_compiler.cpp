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
#include "casadi/core/casadi_logger.hpp"
#include <fstream>
#include <dlfcn.h>
#include <cstdlib>
#include <unistd.h>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_IMPORTER_SHELL_EXPORT
  casadi_register_importer_shell(ImporterInternal::Plugin* plugin) {
    plugin->creator = ShellCompiler::creator;
    plugin->name = "shell";
    plugin->doc = ShellCompiler::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ShellCompiler::options_;
    return 0;
  }

  extern "C"
  void CASADI_IMPORTER_SHELL_EXPORT casadi_load_importer_shell() {
    ImporterInternal::registerPlugin(casadi_register_importer_shell);
  }

  ShellCompiler::ShellCompiler(const std::string& name) :
    ImporterInternal(name) {
      handle_ = 0;
  }

  ShellCompiler::~ShellCompiler() {
    // Unload
    if (handle_) dlclose(handle_);

    if (cleanup_) {
      // Delete the temporary file
      std::string rmcmd = "rm " + bin_name_;
      if (system(rmcmd.c_str())) {
        casadi_warning("Failed to delete temporary file:" + bin_name_);
      }
    }
  }

  Options ShellCompiler::options_
  = {{&ImporterInternal::options_},
     {{"compiler",
       {OT_STRING,
        "Compiler command"}},
      {"linker",
       {OT_STRING,
        "Linker command"}},
      {"folder",
       {OT_STRING,
        "Folder to put temporary objects in."}},
      {"compiler_setup",
       {OT_STRING,
        "Compiler setup command. Intended to be fixed."
        " The 'flag' option is the prefered way to set"
        " custom flags."}},
      {"linker_setup",
       {OT_STRING,
        "Linker setup command. Intended to be fixed."
        " The 'flag' option is the prefered way to set"
        " custom flags."}},
      {"compiler_flags",
       {OT_STRINGVECTOR,
        "Alias for 'compiler_flags'"}},
      {"flags",
        {OT_STRINGVECTOR,
        "Compile flags for the JIT compiler. Default: None"}},
      {"linker_flags",
       {OT_STRINGVECTOR,
        "Linker flags for the JIT compiler. Default: None"}},
      {"cleanup",
       {OT_BOOL,
        "Cleanup temporary files when unloading. Default: true"}}
     }
  };

  void ShellCompiler::init(const Dict& opts) {
    // Base class
    ImporterInternal::init(opts);

    // Default options
    string compiler = "gcc";
    cleanup_ = true;

    string linker = "gcc";
    string compiler_setup = "-fPIC -c";
    string linker_setup = "-shared";
    vector<string> compiler_flags;
    vector<string> linker_flags;


    // Read options
    for (auto&& op : opts) {
      if (op.first=="compiler") {
        compiler = op.second.to_string();
      } else if (op.first=="linker") {
        linker = op.second.to_string();
      } else if (op.first=="compiler_setup") {
        compiler_setup = op.second.to_string();
      } else if (op.first=="cleanup") {
        cleanup_ = op.second;
      } else if (op.first=="linker_setup") {
        linker_setup = op.second.to_string();
      } else if (op.first=="compiler_flags" || op.first=="flags") {
        compiler_flags = op.second;
      } else if (op.first=="linker_flags") {
        linker_flags = op.second;
      }
    }

    // Name of temporary file
#ifdef HAVE_MKSTEMPS
    // Preferred solution
    char obj_name[] = "tmp_casadi_compiler_shell_XXXXXX.o";
    if (mkstemps(obj_name, 2) == -1) {
      casadi_error("Failed to create a temporary object file name");
    }
    obj_name_ = obj_name;

    char bin_name[] = "tmp_casadi_compiler_shell_XXXXXX.so";
    if (mkstemps(bin_name, 3) == -1) {
      casadi_error("Failed to create a temporary library file name");
    }
    bin_name_ = bin_name;
#else
    // Fallback, may result in deprecation warnings
    char* bin_name = tempnam(0, "ca.so");
    bin_name_ = bin_name;
    free(bin_name);

    // Fallback, may result in deprecation warnings
    char* obj_name = tempnam(0, "ca.o");
    obj_name_ = obj_name;
    free(obj_name);
#endif

    // Have relative paths start with ./
    if (obj_name_.at(0)!='/') {
      obj_name_ = "./" + obj_name_;
    }

    if (bin_name_.at(0)!='/') {
      bin_name_ = "./" + bin_name_;
    }

    // Construct the compiler command
    stringstream cccmd;
    cccmd << compiler;
    for (vector<string>::const_iterator i=compiler_flags.begin(); i!=compiler_flags.end(); ++i) {
      cccmd << " " << *i;
    }
    cccmd << " " << compiler_setup;

    // C/C++ source file
    cccmd << " " << name_;

    // Temporary object file
    cccmd << " -o " << obj_name_;

    // Compile into an object
    if (verbose_) userOut() << "calling \"" << cccmd.str() + "\"" << std::endl;
    if (system(cccmd.str().c_str())) {
      casadi_error("Compilation failed. Tried \"" + cccmd.str() + "\"");
    }

    // Link step
    stringstream ldcmd;
    ldcmd << linker;
    for (vector<string>::const_iterator i=linker_flags.begin(); i!=linker_flags.end(); ++i) {
      ldcmd << " " << *i;
    }
    ldcmd << " " << linker_setup;

    // Temporary file
    ldcmd << " " << obj_name_ << " -o " << bin_name_;

    // Compile into a shared library
    if (verbose_) userOut() << "calling \"" << ldcmd.str() << "\"" << std::endl;
    if (system(ldcmd.str().c_str())) {
      casadi_error("Linking failed. Tried \"" + ldcmd.str() + "\"");
    }

    // Load shared library
    handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0,
      "CommonExternal: Cannot open function: " + bin_name_ + ". error code: " +
      str(dlerror()));
    // reset error
    dlerror();
  }

  signal_t ShellCompiler::get_function(const std::string& symname) {
    signal_t ret;
    ret = reinterpret_cast<signal_t>(dlsym(handle_, symname.c_str()));
    if (dlerror()) {
      ret=0;
      dlerror(); // Reset error flags
    }
    return ret;
  }

} // namespace casadi
