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
    if (can_have_meta()) {
      std::vector<std::string> lines;
      int offset;
      get_meta(lines, offset);
      parse(lines, offset);
      if (has_meta("SYMBOLS")) {
        meta_symbols_ = text2set<std::string>(commands_.at("SYMBOLS").second);
      }
    }
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

  bool CompilerInternal::has_function(const std::string& symname) const {
    // Check if in meta information
    if (meta_symbols_.count(symname)) return true;

    // Convert to a dummy function pointer
    return const_cast<CompilerInternal*>(this)->get_function(symname)!=0;
  }

  DllLibrary::DllLibrary(const std::string& bin_name)
    : CompilerInternal(bin_name), handle_(0) {
#ifdef WITH_DL
#ifdef _WIN32
    handle_ = LoadLibrary(TEXT(name_.c_str()));
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open \""
                          << name_ << "\". Error code (WIN32): "<< GetLastError());
#else // _WIN32
    handle_ = dlopen(name_.c_str(), RTLD_LAZY);
    casadi_assert_message(handle_!=0, "CommonExternal: Cannot open \""
                          << name_ << "\". Error code: "<< dlerror());
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

  signal_t DllLibrary::get_function(const std::string& sym) {
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

  std::string CompilerInternal::get_meta(const std::string& cmd, int ind) const {
    if (ind>=0) return get_meta(indexed(cmd, ind));
    casadi_assert_message(has_meta(cmd), "No such command: " + cmd);
    return commands_.at(cmd).second;
  }

  bool CompilerInternal::has_meta(const std::string& cmd, int ind) const {
    if (ind>=0) return has_meta(indexed(cmd, ind));
    return commands_.find(cmd) != commands_.end();
  }

  void CompilerInternal::parse(const std::vector<std::string>& lines, int offset) {
    // Loop over the lines
    auto line_it = lines.cbegin();
    while (line_it!=lines.cend()) {
      // If comment or empty line, skip
      if (line_it->empty() || line_it->at(0)=='#') {
        line_it++;
        continue;
      }

      // Current line number
      int line_no = line_it - lines.cbegin() + 1 + offset;

      // Get command string
      casadi_assert_message(line_it->at(0)==':',
                            "Syntax error: " + *line_it + " is not a command string");
      string cmd = line_it->substr(1, line_it->find(' ')-1);

      // New entry
      stringstream ss;

      // Collect the meta data
      size_t start = cmd.size()+2;
      while (true) {
        // Find the backslash, if any
        size_t stop = line_it->find('\\');

        // Add to entry
        ss << line_it->substr(start, stop-start);

        // Break if no more lines or not multiline
        if (++line_it==lines.cend() || stop == string::npos) break;

        // Multiline entry
        if (start!=stop) ss << std::endl;

        // No offset for subsequent lines
        start = 0;
      }

      // Insert new element in map
      auto new_el = commands_.insert(make_pair(cmd, make_pair(line_no, ss.str())));
      casadi_assert_message(new_el.second, "Duplicate entry: \"" + cmd + "\"");
    }
  }

} // namespace casadi
