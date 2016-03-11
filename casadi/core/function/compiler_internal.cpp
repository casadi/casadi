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


#include "compiler_internal.hpp"

using namespace std;
namespace casadi {

  CompilerInternal::CompilerInternal(const std::string& name) : name_(name) {
  }

  CompilerInternal::~CompilerInternal() {
  }

  void CompilerInternal::print(ostream &stream) const {
    stream << "Compiler" << endl;
  }

  void CompilerInternal::repr(ostream &stream) const {
    stream << "Compiler" << endl;
  }

  std::map<std::string, CompilerInternal::Plugin> CompilerInternal::solvers_;

  const std::string CompilerInternal::infix_ = "compiler";

  Options CompilerInternal::options_
  = {{},
     {{}
     }
  };

  void CompilerInternal::construct(const Dict& opts) {
    // Get a reference to the options structure
    const Options& options = get_options();

    // Make sure all options exist and have the correct type
    for (auto&& op : opts) {
      const Options::Entry* entry = options.find(op.first);
      casadi_assert_message(entry!=0, "No such option: " + op.first);
      casadi_assert_message(op.second.can_cast_to(entry->type),
                            "Illegal type for " + op.first);
    }

    init(opts);
  }

  void CompilerInternal::init(const Dict& opts) {
    // Read meta information from file
    std::vector<std::string> lines;
    int offset;
    get_meta(lines, offset);
    meta_ = ParsedFile(lines, offset);
  }

  void CompilerInternal::get_meta(std::vector<std::string>& lines, int& offset) const {
    // Open source file and search for meta information
    ifstream file(name_);
    std::string line;
    offset = 0;
    lines.clear();
    while (getline(file, line)) {
      // Update offset
      offset++;

      // Try to find a /*CASADI delimiter
      if (line.find("/*CASADIMETA") != string::npos) {
        // Delimimiter found, find */ delimiter
        while (getline(file, line)) {
          if (line.find("*/") != string::npos) return; // Successful return
          lines.push_back(line);
        }
        casadi_error("End-of-file reached while searching for \"*/\"");
      }
    }
  }

  DllLibrary::DllLibrary(const std::string& bin_name)
    : bin_name_(bin_name), handle_(0) {
#ifdef WITH_DL
#ifdef _WIN32
    handle_ = LoadLibrary(TEXT(bin_name_.c_str()));
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open \""
                          << bin_name_ << "\". Error code (WIN32): "<< GetLastError());
#else // _WIN32
    handle_ = dlopen(bin_name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open \""
                          << bin_name_ << "\". Error code: "<< dlerror());
    // reset error
    dlerror();
#endif // _WIN32
#else // WITH_DL
    casadi_error("CommonExternal: WITH_DL  not activated");
#endif // WITH_DL
  }

  DllLibrary::~DllLibrary() {
#ifdef WITH_DL
    // close the dll
#ifdef _WIN32
    if (handle_) FreeLibrary(handle_);
#else // _WIN32
    if (handle_) dlclose(handle_);
#endif // _WIN32
#endif // WITH_DL
  }

  bool DllLibrary::has(const std::string& sym) {
    return get(sym)!=0;
  }

  signal_t DllLibrary::get(const std::string& sym) {
#ifdef WITH_DL
#ifdef _WIN32
    return (signal_t)GetProcAddress(handle_, TEXT(sym.c_str()));
#else // _WIN32
    signal_t fcnPtr = (signal_t)dlsym(handle_, sym.c_str());
    if (dlerror()) {
      fcnPtr=0;
      dlerror(); // Reset error flags
    }
    return fcnPtr;
#endif // _WIN32
#endif // WITH_DL
  }

  const ParsedFile& DllLibrary::meta() const {
    static ParsedFile singleton;
    return singleton;
  }

  JitLibrary::JitLibrary(const Compiler& compiler)
    : compiler_(compiler) {
    if (compiler.meta().has("SYMBOLS")) {
      meta_symbols_ = compiler.meta().to_set<std::string>("SYMBOLS");
    }
  }

  JitLibrary::~JitLibrary() {
  }

  bool JitLibrary::has(const std::string& sym) {
    // Check if in meta information
    if (meta_symbols_.count(sym)) return true;

    // Convert to a dummy function pointer
    return get(sym)!=0;
  }

  signal_t JitLibrary::get(const std::string& sym) {
    return (signal_t)compiler_.get_function(sym);
  }

  const ParsedFile& JitLibrary::meta() const {
    return compiler_.meta();
  }

} // namespace casadi
