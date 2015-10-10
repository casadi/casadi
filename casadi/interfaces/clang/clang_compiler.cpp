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

    myerr_ = 0;
    executionEngine_ = 0;
    context_ = 0;
    act_ = 0;
  }

  ClangCompiler::~ClangCompiler() {
    if (act_) delete act_;
    if (myerr_) delete myerr_;
    if (executionEngine_) delete executionEngine_;
    if (context_) delete context_;
  }

  void ClangCompiler::init() {
    // Initialize the base classes
    CompilerInternal::init();

    // Arguments to pass to the clang frontend
    vector<const char *> args(1, name_.c_str());
    std::vector<std::string> flags;
    if (hasSetOption("flags")) {
      flags = getOption("flags");
      for (auto i=flags.begin(); i!=flags.end(); ++i) {
        args.push_back(i->c_str());
      }
    }

    // Create the compiler instance
    clang::CompilerInstance compInst;

    // A symbol in the DLL
    void *addr = reinterpret_cast<void*>(&casadi_register_compiler_clang);

    // Get runtime include path
    std::string jit_include, filesep;
#ifdef _WIN32
    char buffer[MAX_PATH];
    HMODULE hm = NULL;
    if (!GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            (LPCSTR)addr, &hm)) {
      casadi_error("GetModuleHandle failed");
    }
    GetModuleFileNameA(hm, buffer, sizeof(buffer));
    PathRemoveFileSpecA(buffer);
    jit_include = buffer;
    filesep = "\\";
#else // _WIN32
    Dl_info dl_info;
    if (!dladdr(addr, &dl_info)) {
      casadi_error("dladdr failed");
    }
    jit_include = dl_info.dli_fname;
    jit_include = jit_include.substr(0, jit_include.find_last_of('/'));
    filesep = "/";
#endif // _WIN32
    jit_include += filesep + "casadi" + filesep + "jit";

#if 0
    // Initialize target info with the default triple for our platform.
    auto targetoptions = std::make_shared<clang::TargetOptions>();
    targetoptions->Triple = llvm::sys::getDefaultTargetTriple();
    clang::TargetInfo *targetInfo =
      clang::TargetInfo::CreateTargetInfo(compInst.getDiagnostics(), targetoptions);
    compInst.setTarget(targetInfo);
#endif

    // The compiler invocation needs a DiagnosticsEngine so it can report problems
    clang::DiagnosticOptions* diagOpts = new clang::DiagnosticOptions();
    myerr_ = new llvm::raw_os_ostream(userOut<true>());
    clang::TextDiagnosticPrinter *diagClient = new clang::TextDiagnosticPrinter(*myerr_, diagOpts);

    clang::DiagnosticIDs* diagID = new clang::DiagnosticIDs();
    // This object takes ownerships of all three passed-in pointers
    clang::DiagnosticsEngine diags(diagID, diagOpts, diagClient);

    // Create the compiler invocation
    clang::CompilerInvocation* compInv = new clang::CompilerInvocation();
    clang::CompilerInvocation::CreateFromArgs(*compInv, &args[0],
                                              &args[0] + args.size(), diags);
    compInst.setInvocation(compInv);

    // Get ready to report problems
    compInst.createDiagnostics();
    if (!compInst.hasDiagnostics())
      casadi_error("Cannot create diagnostics");

    // Set resource directory
    std::string resourcedir = jit_include + filesep + "clang" + filesep + CLANG_VERSION_STRING;
    compInst.getHeaderSearchOpts().ResourceDir = resourcedir;

    // Read the system includes (C or C++)
    vector<pair<string, bool> > system_include = getIncludes("system_includes.txt", jit_include);
    for (auto i=system_include.begin(); i!=system_include.end(); ++i) {
      compInst.getHeaderSearchOpts().AddPath(i->first,
                                             clang::frontend::System, i->second, false);
    }

    // Read the system includes (C only)
    system_include = getIncludes("csystem_includes.txt", jit_include);
    for (auto i=system_include.begin(); i!=system_include.end(); ++i) {
      compInst.getHeaderSearchOpts().AddPath(i->first,
                                             clang::frontend::CSystem, i->second, false);
    }

    // Read the system includes (C++ only)
    system_include = getIncludes("cxxsystem_includes.txt", jit_include);
    for (auto i=system_include.begin(); i!=system_include.end(); ++i) {
      compInst.getHeaderSearchOpts().AddPath(i->first,
                                             clang::frontend::CXXSystem, i->second, false);
    }

#ifdef _WIN32
    char pathsep = ';';
#else
    char pathsep = ':';
#endif

    // Search path
    std::stringstream paths;
    paths << getOption("include_path").toString() << pathsep;
    std::string path;
    while (std::getline(paths, path, pathsep)) {
      compInst.getHeaderSearchOpts().AddPath(path.c_str(), clang::frontend::System, false, false);
    }

    // Create an LLVM context (NOTE: should use a static context instead?)
    context_ = new llvm::LLVMContext();

    // Create an action and make the compiler instance carry it out
    act_ = new clang::EmitLLVMOnlyAction(context_);
    if (!compInst.ExecuteAction(*act_))
      casadi_error("Cannot execute action");

    // Grab the module built by the EmitLLVMOnlyAction
    #if LLVM_VERSION_MAJOR>=3 && LLVM_VERSION_MINOR>=5
    std::unique_ptr<llvm::Module> module = act_->takeModule();
    module_ = module.get();
    #else
    llvm::Module* module = act_->takeModule();
    module_ = module;
    #endif

    llvm::InitializeNativeTarget();
    llvm::InitializeNativeTargetAsmPrinter();

    // Create the JIT.  This takes ownership of the module.
    std::string ErrStr;
    executionEngine_ =
      llvm::EngineBuilder(std::move(module)).setEngineKind(llvm::EngineKind::JIT)
      .setErrorStr(&ErrStr).create();
    if (!executionEngine_) {
      casadi_error("Could not create ExecutionEngine: " << ErrStr);
    }

    executionEngine_->finalizeObject();
  }

  void* ClangCompiler::getFunction(const std::string& symname) {
    return reinterpret_cast<void*>((intptr_t)executionEngine_
                                   ->getPointerToFunction(module_->getFunction(symname)));
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
    
    userOut() << "test001c" << std::endl;
    std::string line;
    userOut() << "test002" << std::endl;

    while (std::getline(setup_file, line)) {
      // Skip empty lines
      if (line.empty()) continue;

      userOut() << "test003" << std::endl;

      // Check if framework
      size_t loc = line.find(" (framework directory)");
      bool isframework = loc != string::npos;

      userOut() << "test004" << std::endl;

      if (isframework) {
        // Truncate path
        line = line.substr(0, loc);
      }

      userOut() << "test005" << std::endl;

      // Check if the path is absolute or relative
#ifdef _WIN32
      bool relative = PathIsRelative(TEXT(line.c_str()));
#else // _WIN32
      bool relative = line.at(0)!=sep;
#endif // _WIN32

      userOut() << "test006" << std::endl;


      if (relative) {
        // Relative path, make absolute
        ret.push_back(make_pair(path + sep + line, isframework));
      } else {
        // Absolute path
        ret.push_back(make_pair(line, isframework));
      }

      userOut() << "test007" << std::endl;

    }

    return ret;
  }

} // namespace casadi
