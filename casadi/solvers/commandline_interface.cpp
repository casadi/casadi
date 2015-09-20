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


#include "commandline_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/casadi_meta.hpp"
#include <fstream>
#include <stdlib.h>
#include <dlfcn.h>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_JITCOMPILER_COMMANDLINE_EXPORT
  casadi_register_jitcompiler_commandline(JitCompilerInternal::Plugin* plugin) {
    plugin->creator = CommandlineJitCompiler::creator;
    plugin->name = "commandline";
    plugin->doc = CommandlineJitCompiler::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_JITCOMPILER_COMMANDLINE_EXPORT casadi_load_jitcompiler_commandline() {
    JitCompilerInternal::registerPlugin(casadi_register_jitcompiler_commandline);
  }

  CommandlineJitCompiler* CommandlineJitCompiler::clone() const {
    // Return a deep copy
    CommandlineJitCompiler* node = new CommandlineJitCompiler(name_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  CommandlineJitCompiler::CommandlineJitCompiler(const std::string& name) :
    JitCompilerInternal(name) {
    addOption("compiler", OT_STRING, "gcc", "Compiler command");
    addOption("compiler_setup", OT_STRING, "-fPIC -shared", "Compiler setup command");
    addOption("flags", OT_STRINGVECTOR, GenericType(),
      "Compile flags for the JIT compiler. Default: None");
  }

  CommandlineJitCompiler::~CommandlineJitCompiler() {
    // Unload
    if (handle_) dlclose(handle_);

    // Delete the temporary file
    std::string cmd = "rm " + bin_name_;
    int flag = system(cmd.c_str());
    casadi_assert_warning(flag==0, "Failed to delete temporary file");
  }

  void CommandlineJitCompiler::init() {
    // Initialize the base classes
    JitCompilerInternal::init();

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
    char bin_name[] = "tmp_casadi_jitcompiler_commandline_XXXXXX.so";
    int flag = mkstemps(bin_name, 3);
    bin_name_ = bin_name;
    cmd << " -o " << bin_name_;

    // Compile into a shared library
    flag = system(cmd.str().c_str());
    if (flag!=0) {
      casadi_error("Compilation failed. Tried \"" + cmd.str() + "\"");
    }

    // Load shared library
    handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open function: "
                          << bin_name_ << ". error code: "<< dlerror());
    // reset error
    dlerror();
  }

  void* CommandlineJitCompiler::getFunction(const std::string& symname) {
    void* ret;
    ret = reinterpret_cast<void*>(dlsym(handle_, symname.c_str()));
    if (dlerror()) {
      ret=0;
      dlerror(); // Reset error flags
    }
    return ret;
  }

} // namespace casadi
