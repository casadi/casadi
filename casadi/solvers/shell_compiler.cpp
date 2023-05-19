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


#include "shell_compiler.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_meta.hpp"
#include "casadi/core/casadi_logger.hpp"
#include <fstream>

// Set default object file suffix
#ifndef OBJECT_FILE_SUFFIX
#define OBJECT_FILE_SUFFIX CASADI_OBJECT_FILE_SUFFIX
#endif // OBJECT_FILE_SUFFIX

#include <cstdlib>

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
      handle_ = nullptr;
  }

  ShellCompiler::~ShellCompiler() {
    if (handle_) close_shared_library(handle_);

    if (cleanup_) {
      if (remove(bin_name_.c_str())) casadi_warning("Failed to remove " + bin_name_);
      if (remove(obj_name_.c_str())) casadi_warning("Failed to remove " + obj_name_);
      for (const std::string& s : extra_suffixes_) {
        std::string name = base_name_+s;
        remove(name.c_str());
      }
    }
  }

  const Options ShellCompiler::options_
  = {{&ImporterInternal::options_},
     {{"compiler",
       {OT_STRING,
        "Compiler command"}},
      {"linker",
       {OT_STRING,
        "Linker command"}},
      {"directory",
       {OT_STRING,
        "Directory to put temporary objects in. Must end with a file separator."}},
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
        "Cleanup temporary files when unloading. Default: true"}},
      {"compiler_output_flag",
       {OT_STRING,
       "Compiler flag to denote object output. Default: '-o '"}},
      {"linker_output_flag",
       {OT_STRING,
       "Linker flag to denote shared library output. Default: '-o '"}},
      {"extra_suffixes",
       {OT_STRINGVECTOR,
       "List of suffixes for extra files that the compiler may generate. Default: None"}},
      {"name",
       {OT_STRING,
        "The file name used to write out compiled objects/libraries. "
        "The actual file names used depend on 'temp_suffix' and include extensions. "
        "Default: 'tmp_casadi_compiler_shell'"}},
      {"temp_suffix",
       {OT_BOOL,
        "Use a temporary (seemingly random) filename suffix for file names. "
        "This is desired for thread-safety. "
        "This behaviour may defeat caching compiler wrappers. "
        "Default: true"}},
     }
  };

  void ShellCompiler::init(const Dict& opts) {
    // Base class
    ImporterInternal::init(opts);

    // Default options

    cleanup_ = true;
    bool temp_suffix = true;
    std::string bare_name = "tmp_casadi_compiler_shell";
    std::string directory = "";


    std::vector<std::string> compiler_flags;
    std::vector<std::string> linker_flags;
    std::string suffix = OBJECT_FILE_SUFFIX;

    std::string compiler = GlobalOptions::default_compiler;
    std::string linker = GlobalOptions::default_linker;
    std::string compiler_setup = GlobalOptions::default_compiler_setup;
    std::string linker_setup = GlobalOptions::default_linker_setup;
    std::string compiler_output_flag = GlobalOptions::default_compiler_output_flag;
    std::string linker_output_flag = GlobalOptions::default_linker_output_flag;
    extra_suffixes_ = GlobalOptions::default_compiler_extra_suffixes;
  
    // Read options
    for (auto&& op : opts) {
      if (op.first=="compiler") {
        compiler = op.second.to_string();
      } else if (op.first=="linker") {
        linker = op.second.to_string();
      } else if (op.first=="directory") {
        directory = op.second.to_string();
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
      } else if (op.first=="compiler_output_flag") {
        compiler_output_flag = op.second.to_string();
      } else if (op.first=="linker_output_flag") {
        linker_output_flag = op.second.to_string();
      } else if (op.first=="extra_suffixes") {
        extra_suffixes_ = op.second.to_string_vector();
      } else if (op.first=="name") {
        bare_name = op.second.to_string();
      } else if (op.first=="temp_suffix") {
        temp_suffix = op.second;
      }
    }

    // Name of temporary file
    if (temp_suffix) {
      obj_name_ = temporary_file(directory + bare_name, suffix);
    } else {
      obj_name_ = directory + bare_name + suffix;
    }
    base_name_ = std::string(obj_name_.begin(), obj_name_.begin()+obj_name_.size()-suffix.size());
    bin_name_ = base_name_+SHARED_LIBRARY_SUFFIX;

#ifndef _WIN32
    // Have relative paths start with ./
    if (obj_name_.at(0)!='/') {
      obj_name_ = "./" + obj_name_;
    }

    if (bin_name_.at(0)!='/') {
      bin_name_ = "./" + bin_name_;
    }
#endif // _WIN32

    // Construct the compiler command
    std::stringstream cccmd;
    cccmd << compiler;
    for (auto i=compiler_flags.begin(); i!=compiler_flags.end(); ++i) {
      cccmd << " " << *i;
    }
    cccmd << " " << compiler_setup;

    // C/C++ source file
    cccmd << " " << name_;

    // Temporary object file
    cccmd << " " + compiler_output_flag << obj_name_;

    // Compile into an object
    if (verbose_) casadi_message("calling \"" + cccmd.str() + "\"");
    if (system(cccmd.str().c_str())) {
      casadi_error("Compilation failed. Tried \"" + cccmd.str() + "\"");
    }

    // Link step
    std::stringstream ldcmd;
    ldcmd << linker;

    // Temporary file
    ldcmd << " " << obj_name_ << " " + linker_output_flag + bin_name_;

    // Add flags
    for (auto i=linker_flags.begin(); i!=linker_flags.end(); ++i) {
      ldcmd << " " << *i;
    }
    ldcmd << " " << linker_setup;

    // Compile into a shared library
    if (verbose_) casadi_message("calling \"" + ldcmd.str() + "\"");
    if (system(ldcmd.str().c_str())) {
      casadi_error("Linking failed. Tried \"" + ldcmd.str() + "\"");
    }

    std::vector<std::string> search_paths = get_search_paths();
    handle_ = open_shared_library(bin_name_, search_paths, "ShellCompiler::init");

  }

  std::string ShellCompiler::library() const {
    return bin_name_;
  }

  signal_t ShellCompiler::get_function(const std::string& symname) {
#ifdef _WIN32
    return (signal_t)GetProcAddress(handle_, TEXT(symname.c_str()));
#else // _WIN32
    signal_t fcnPtr = reinterpret_cast<signal_t>(dlsym(handle_, symname.c_str()));
    if (dlerror()) {
      fcnPtr=nullptr;
      dlerror(); // Reset error flags
    }
    return fcnPtr;
#endif // _WIN32
  }

} // namespace casadi
