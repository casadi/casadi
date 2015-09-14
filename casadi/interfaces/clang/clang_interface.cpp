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


#include "clang_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/function/mx_function.hpp"
#include "casadi/core/casadi_meta.hpp"

using namespace clang;
using namespace clang::driver;


using namespace std;
namespace casadi {

  extern "C"
  int CASADI_JITCOMPILER_CLANG_EXPORT
  casadi_register_jitcompiler_clang(JitCompilerInternal::Plugin* plugin) {
    plugin->creator = ClangJitCompilerInterface::creator;
    plugin->name = "clang";
    plugin->doc = ClangJitCompilerInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_JITCOMPILER_CLANG_EXPORT casadi_load_jitcompiler_clang() {
    JitCompilerInternal::registerPlugin(casadi_register_jitcompiler_clang);
  }

  ClangJitCompilerInterface* ClangJitCompilerInterface::clone() const {
    // Return a deep copy
    ClangJitCompilerInterface* node = new ClangJitCompilerInterface(f_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  ClangJitCompilerInterface::ClangJitCompilerInterface(const Function& f) :
      JitCompilerInternal(f) {
    addOption("include_path", OT_STRING, "", "Include paths for the JIT compiler."
      " The include directory shipped with CasADi will be automatically appended.");
    addOption("flags", OT_STRINGVECTOR, GenericType(),
      "Compile flags for the JIT compiler. Default: -O3.");
  }

  ClangJitCompilerInterface::~ClangJitCompilerInterface() {
    if (isInit()) {
      if (TheExecutionEngine) {
        delete TheExecutionEngine;
        TheExecutionEngine=0;
      }
      if (Context) {
        delete Context;
        Context=0;
      }
    }
    //llvm::llvm_shutdown();

  }

  void ClangJitCompilerInterface::init() {
    // Initialize the base classes
    JitCompilerInternal::init();
    log("ClangJitCompilerInterface::init", "Enter");



    std::string name = "foo";

    f_.generate(name);

    // Path to the C file
    string inputPath = name + ".c";


    Context = new llvm::LLVMContext();

    // Arguments to pass to the clang frontend
    vector<const char *> args;
    args.push_back(inputPath.c_str());

    std::vector<std::string> flags;

    if (hasSetOption("flags")) {
      flags = getOption("flags");
      for (int i=0;i<flags.size();++i) {
        args.push_back(flags[i].c_str());
      }
    } else {
      args.push_back("-O3");
    }

    // The compiler invocation needs a DiagnosticsEngine so it can report problems
    DiagnosticOptions* DiagOpts = new DiagnosticOptions();
    llvm::raw_ostream* myerr = new llvm::raw_os_ostream(userOut<true>());
    TextDiagnosticPrinter *DiagClient = new TextDiagnosticPrinter(*myerr, DiagOpts);

    DiagnosticIDs* DiagID = new DiagnosticIDs();
    // This object takes ownerships of all three passed-in pointers
    DiagnosticsEngine Diags(DiagID, DiagOpts, DiagClient);

    // Create the compiler invocation
    clang::CompilerInvocation* CI = new clang::CompilerInvocation();
    clang::CompilerInvocation::CreateFromArgs(*CI, &args[0], &args[0] + args.size(), Diags);

    vector<string> argsFromInvocation;

    // Create the compiler instance
    clang::CompilerInstance Clang;
    Clang.setInvocation(CI);

    // Get ready to report problems
    Clang.createDiagnostics();
    if (!Clang.hasDiagnostics())
      casadi_error("Cannot create diagnostics");

    #ifdef _WIN32
    char pathsep = ';';
    const std::string filesep("\\");
    #else
    char pathsep = ':';
    const std::string filesep("/");
    #endif

    // Search path
    std::stringstream paths;
    paths << getOption("include_path").toString() << pathsep;
    paths << CasadiOptions::getJitIncludePath() << pathsep;
    paths << CasadiMeta::getInstallPrefix() << filesep <<
      "include" << filesep << "casadi" << filesep << "jit";
    std::string path;
    while (std::getline(paths, path, pathsep)) {
      log("ClangInterface::init", std::string("Adding path:") + path);
      Clang.getHeaderSearchOpts().AddPath(path.c_str(), frontend::CSystem, false, false);
    }

    // Create an action and make the compiler instance carry it out
    Act = new clang::EmitLLVMOnlyAction(Context);
    if (!Clang.ExecuteAction(*Act))
      casadi_error("Cannot execute action");

    // Grab the module built by the EmitLLVMOnlyAction
    std::unique_ptr<llvm::Module> module = Act->takeModule();

    llvm::InitializeNativeTarget();
    llvm::InitializeNativeTargetAsmPrinter();
    llvm::Module &M = *module;

    // Create the JIT.  This takes ownership of the module.
    std::string ErrStr;
    TheExecutionEngine =
      llvm::EngineBuilder(std::move(module)).setEngineKind(llvm::EngineKind::JIT)
        .setErrorStr(&ErrStr).create();
    if (!TheExecutionEngine) {
      casadi_error("Could not create ExecutionEngine: " << ErrStr);
    }

    TheExecutionEngine->finalizeObject();

    delete Act; Act = 0;


    std::string name_init = name + "_init";
    std::string name_sparsity = name + "_sparsity";
    std::string name_work = name + "_work";
    std::string name_eval = name;

    init_fun_ = initPtr(intptr_t(
      TheExecutionEngine->getPointerToFunction(M.getFunction(name_init))));
    if (!init_fun_) casadi_error("Symbol " << name_init << "not found.");

    int f_type, n_in, n_out, sz_arg, sz_res;
    int myres = init_fun_(&f_type, &n_in, &n_out, &sz_arg, &sz_res);
    casadi_assert_message(myres==0, "ClangJitCompilerInternal: \"init\" failed");

    ibuf_.resize(n_in);
    obuf_.resize(n_out);

    sparsity_fun_ = sparsityPtr(intptr_t(
      TheExecutionEngine->getPointerToFunction(M.getFunction(name_sparsity))));
    if (!sparsity_fun_) casadi_error("Symbol " << name_sparsity << "not found.");
    work_fun_ = workPtr(intptr_t(
      TheExecutionEngine->getPointerToFunction(M.getFunction(name_work))));
    if (!work_fun_) casadi_error("Symbol " << name_work << "not found.");
    eval_fun_ = evalPtr(intptr_t(
      TheExecutionEngine->getPointerToFunction(M.getFunction(name_eval))));
    if (!eval_fun_) casadi_error("Symbol " << name_eval << "not found.");

    // Get the sparsity patterns
    for (int i=0; i<nIn()+nOut(); ++i) {
      // Get sparsity from file
      int nrow, ncol;
      const int *colind, *row;
      int flag = sparsity_fun_(i, &nrow, &ncol, &colind, &row);
      casadi_assert_message(flag==0, "ClangJitCompilerInternal: \"sparsity\" failed");

      // Col offsets
      vector<int> colindv(colind, colind+ncol+1);

      // Number of nonzeros
      int nnz = colindv.back();

      // Rows
      vector<int> rowv(row, row+nnz);

      // Sparsity
      Sparsity sp(nrow, ncol, colindv, rowv);

      // Save to inputs/outputs
      if (i<nIn()) {
        input(i) = Matrix<double>::zeros(sp);
      } else {
        output(i-nIn()) = Matrix<double>::zeros(sp);
      }
    }

    delete myerr;myerr=0;

    int n_iw, n_w;
    int flag = work_fun_(&n_iw, &n_w);
    casadi_assert_message(flag==0, "ClangJitCompilerInternal: \"work\" failed");
    alloc_iw(n_iw);
    alloc_w(n_w);

    alloc_arg(max(n_in, sz_arg));
    alloc_res(max(n_out, sz_res));
  }

  void ClangJitCompilerInterface::evalD(const double** arg, double** res, int* iw, double* w) {
    eval_fun_(arg, res, iw, w);
  }


} // namespace casadi
