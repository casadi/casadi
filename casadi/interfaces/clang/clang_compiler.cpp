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


#include "clang_compiler.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/casadi_meta.hpp"
#include <fstream>

// To be able to get the plugin path
#ifdef _WIN32 // also for 64-bit
#define NOMINMAX
#include <windows.h>
#include <shlwapi.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_COMPILER_CLANG_EXPORT
  casadi_register_compiler_clang(CompilerInternal::Plugin* plugin) {
    plugin->creator = ClangCompiler::creator;
    plugin->name = "clang";
    plugin->doc = ClangCompiler::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_COMPILER_CLANG_EXPORT casadi_load_compiler_clang() {
    CompilerInternal::registerPlugin(casadi_register_compiler_clang);
  }

  ClangCompiler* ClangCompiler::clone() const {
    // Return a deep copy
    ClangCompiler* node = new ClangCompiler(name_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  ClangCompiler::ClangCompiler(const std::string& name) :
    CompilerInternal(name) {
    addOption("include_path", OT_STRING, "", "Include paths for the JIT compiler."
      " The include directory shipped with CasADi will be automatically appended.");
    addOption("flags", OT_STRINGVECTOR, GenericType(),
      "Compile flags for the JIT compiler. Default: None");

  }

  ClangCompiler::~ClangCompiler() {

  }

  void ClangCompiler::init() {
    // Initialize the base classes
    CompilerInternal::init();
    
    getIncludes("system_includes.txt", "E:\\casadi-matlabR2014b-cad02fe\\casadi\\jit");

  }

  void* ClangCompiler::getFunction(const std::string& symname) {
    
  }

  std::vector<std::pair<std::string, bool> > ClangCompiler::
  getIncludes(const std::string& file, const std::string& path) {
    // File separator
#ifdef _WIN32
    const char sep = '\\';
#else // _WIN32
    const char sep = '/';
#endif // _WIN32

    userOut() << "test001" << std::endl;

    // Return value
    vector<pair<string, bool> > ret;
    userOut() << "test001b" << std::endl;
    
    // Read line-by-line
    std::string file_name = path + sep + file;
    userOut() << "test001cc" << file_name << std::endl;
    std::ifstream setup_file(file_name.c_str());

    userOut() << setup_file.is_open() << std::endl;
    userOut() << setup_file.good() << std::endl;
    userOut() << setup_file.eof() << std::endl;
    userOut() << setup_file.fail() << std::endl;
    userOut() << setup_file.bad() << std::endl;
    userOut() << setup_file.rdstate() << std::endl;
    return ret;
  }

} // namespace casadi
