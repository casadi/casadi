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
#define OBJECT_FILE_SUFFIX CasadiMeta::object_file_suffix()
#endif // OBJECT_FILE_SUFFIX

#include <cstdlib>
#include <cstring>

#ifdef _WIN32
#include <windows.h>
#endif // _WIN32

namespace casadi {

namespace {
#ifdef _WIN32
    // Run a command line; -1 if it could not be started
    int win_run(const std::string& app, const std::string& cmd) {
      STARTUPINFOA si;
      ZeroMemory(&si, sizeof(si));
      si.cb = sizeof(si);
      // No console window when the parent has none
      DWORD flags = GetConsoleCP() ? 0 : CREATE_NO_WINDOW;
      std::vector<char> buf(cmd.begin(), cmd.end());
      buf.push_back('\0');
      PROCESS_INFORMATION pi;
      if (!CreateProcessA(app.c_str(), buf.data(), nullptr, nullptr, TRUE, flags,
          nullptr, nullptr, &si, &pi)) return -1;
      WaitForSingleObject(pi.hProcess, INFINITE);
      DWORD ret = 1;
      GetExitCodeProcess(pi.hProcess, &ret);
      CloseHandle(pi.hProcess);
      CloseHandle(pi.hThread);
      return static_cast<int>(ret);
    }

    // Command interpreter of system()
    std::string win_shell() {
      const char* comspec = getenv("ComSpec");
      if (comspec && *comspec) return comspec;
      char dir[MAX_PATH];
      UINT n = GetSystemDirectoryA(dir, MAX_PATH);
      return std::string(dir, n<MAX_PATH ? n : 0) + "\\cmd.exe";
    }

    // Text as cmd.exe sees it: decoded as by CreateProcessA, optionally with %VAR% expanded
    std::wstring win_wide(const std::string& s, bool expand) {
      std::wstring w(s.size()+1, L'\0');
      int n = MultiByteToWideChar(CP_ACP, 0, s.c_str(), -1, &w[0], static_cast<int>(w.size()));
      w.resize(n>0 ? n-1 : 0);
      if (!expand) return w;
      DWORD m = ExpandEnvironmentStringsW(w.c_str(), nullptr, 0);
      std::wstring e(m, L'\0');
      if (!m || ExpandEnvironmentStringsW(w.c_str(), &e[0], m)!=m) return w;
      e.resize(m-1);
      return e;
    }

    // Back to the code page of CreateProcessA; false if not representable
    bool win_narrow(const std::wstring& w, std::string& s) {
      BOOL lossy = FALSE;
      bool utf8 = GetACP()==CP_UTF8;
      int n = WideCharToMultiByte(CP_ACP, utf8 ? 0 : WC_NO_BEST_FIT_CHARS, w.c_str(), -1,
        nullptr, 0, nullptr, utf8 ? nullptr : &lossy);
      if (n<=0 || lossy) return false;
      s.resize(n);
      WideCharToMultiByte(CP_ACP, utf8 ? 0 : WC_NO_BEST_FIT_CHARS, w.c_str(), -1,
        &s[0], n, nullptr, nullptr);
      s.resize(n-1);
      return true;
    }

    // Does cmd.exe end up outside quotes after an option?
    bool cmd_balanced(const std::wstring& s) {
      bool in_quotes = false;
      for (size_t i=0; i<s.size(); ++i) {
        if (s[i]==L'"') {
          in_quotes = !in_quotes;
        } else if (!in_quotes && s[i]==L'^' && ++i==s.size()) {
          return false;
        }
      }
      return !in_quotes;
    }

    // Do non-ASCII bytes reach a narrow program unchanged, without best fit or trail bytes?
    bool win_bytes_safe() {
      // The program uses the code page of the system, which may differ from this process
      char acp[16] = "";
      DWORD size = sizeof(acp);
      HKEY key;
      if (RegOpenKeyExA(HKEY_LOCAL_MACHINE, "SYSTEM\\CurrentControlSet\\Control\\Nls\\CodePage",
          0, KEY_QUERY_VALUE, &key)!=ERROR_SUCCESS) return false;
      LONG ret = RegQueryValueExA(key, "ACP", nullptr, nullptr,
        reinterpret_cast<BYTE*>(acp), &size);
      RegCloseKey(key);
      if (ret!=ERROR_SUCCESS || size>=sizeof(acp)) return false;
      if (GetACP()==CP_UTF8) return std::string(acp)=="65001";
      CPINFO cp;
      return std::string(acp)==std::to_string(GetACP()) && GetCPInfo(CP_ACP, &cp)
        && cp.MaxCharSize==1;
    }
#endif // _WIN32

    // Quote a path for the shell
    std::string quote(const std::string& s) {
      // Not to be taken for an option
      if (!s.empty() && (s[0]=='-' || s[0]=='@')) return quote("./" + s);
#ifdef _WIN32
      // cmd.exe expands % and !, and with delayed expansion removes ^, also within quotes
      for (char c : s) {
        casadi_assert(c!='%' && c!='!' && c!='^' && c!='"' && static_cast<unsigned char>(c)>=0x20,
          "Unsupported character in '" + s + "'");
      }
      // Non-ASCII could reach the program as " or \, by best fit or as a trail byte
      if (!win_bytes_safe()) {
        for (char c : s) {
          casadi_assert(static_cast<unsigned char>(c)<0x80, "Unsupported character in '" + s + "'");
        }
      }
      // Plain paths as before
      if (!s.empty() && s.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-+/.:@\\~#$")==std::string::npos) return s;
      // Trailing backslashes doubled for the MSVC runtime
      size_t n_backslash = s.size() - (s.find_last_not_of('\\')+1);
      return "\"" + s + std::string(n_backslash, '\\') + "\"";
#else // _WIN32
      // Still interpreted by the shell
      return s;
#endif // _WIN32
    }

    // Include directory, with enclosing quotes removed as by the shell
    std::string include_dir(const std::string& entry) {
      // Without surrounding whitespace, which separated it for the shell
      const char* ws = " \t\r\n";
      size_t first = entry.find_first_not_of(ws), last = entry.find_last_not_of(ws);
      std::string s = first==std::string::npos ? "" : entry.substr(first, last-first+1);
#ifdef _WIN32
      std::string e;
      if (!win_narrow(win_wide(s, true), e)) e = s;
      if (e.size()>=2 && e.front()=='"' && e.back()=='"') return e.substr(1, e.size()-2);
      return e;
#else // _WIN32
      // Still interpreted by the shell
      return s;
#endif // _WIN32
    }

    // Run a command through the shell
    int run(const std::string& cmd) {
#ifdef _WIN32
      // As system(), without a console window
      std::string shell = win_shell();
      return win_run(shell, "\"" + shell + "\" /s /c \"" + cmd + "\"");
#else // _WIN32
      return system(cmd.c_str());
#endif // _WIN32
    }
}  // namespace

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
        " The 'flag' option is the preferred way to set"
        " custom flags."}},
      {"linker_setup",
       {OT_STRING,
        "Linker setup command. Intended to be fixed."
        " The 'flag' option is the preferred way to set"
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
      {"compiler_include_flag",
       {OT_STRING,
       "Compiler flag to add an include directory. Default: '-I' ('/I' on MSVC)"}},
      {"include_dirs",
       {OT_STRINGVECTOR,
       "List of include directories, each passed to the compiler prefixed with "
       "'compiler_include_flag'. OS-agnostic alternative to passing '-I...' via 'flags'. "
       "Default: None"}},
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
    std::string directory = FunctionInternal::get_jit_directory(opts);

    std::vector<std::string> compiler_flags;
    std::vector<std::string> linker_flags;
    std::vector<std::string> include_dirs;
    std::string suffix = OBJECT_FILE_SUFFIX;

#ifdef _WIN32
    std::string compiler = "cl.exe";
    std::string linker = "link.exe";
    std::string compiler_setup = "/c";
    std::string linker_setup = "/DLL";
    std::string compiler_output_flag = "/Fo";
    std::string linker_output_flag = "/out:";
    std::string compiler_include_flag = "/I";
    extra_suffixes_ = {".exp", ".lib"};
#elif defined(__APPLE__)
    std::string compiler = "clang";
    std::string linker = "clang";
    std::string compiler_setup = "-fPIC -c";
    std::string linker_setup = "-shared";
    std::string compiler_output_flag = "-o ";
    std::string linker_output_flag = "-o ";
    std::string compiler_include_flag = "-I";
#else
    std::string compiler = "gcc";
    std::string linker = "gcc";
    std::string compiler_setup = "-fPIC -c";
    std::string linker_setup = "-shared";
    std::string compiler_output_flag = "-o ";
    std::string linker_output_flag = "-o ";
    std::string compiler_include_flag = "-I";
#endif

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
      } else if (op.first=="compiler_output_flag") {
        compiler_output_flag = op.second.to_string();
      } else if (op.first=="linker_output_flag") {
        linker_output_flag = op.second.to_string();
      } else if (op.first=="compiler_include_flag") {
        compiler_include_flag = op.second.to_string();
      } else if (op.first=="include_dirs") {
        include_dirs = op.second.to_string_vector();
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
      obj_name_ = temporary_file(bare_name, suffix, directory);
    } else {
      obj_name_ = directory + bare_name + suffix;
    }
    base_name_ = std::string(obj_name_.begin(), obj_name_.begin()+obj_name_.size()-suffix.size());
    bin_name_ = base_name_+SHARED_LIBRARY_SUFFIX;

#ifdef _WIN32
    std::vector<std::string> all = {compiler, compiler_setup, compiler_output_flag,
      compiler_include_flag, linker, linker_setup, linker_output_flag};
    all.insert(all.end(), compiler_flags.begin(), compiler_flags.end());
    all.insert(all.end(), linker_flags.begin(), linker_flags.end());
    for (const std::string& s : all) {
      for (char c : s) {
        casadi_assert(c=='\t' || static_cast<unsigned char>(c)>=0x20,
          "Unsupported character in '" + s + "'");
      }
      // CasADi's paths must not end up between quotes
      casadi_assert(cmd_balanced(win_wide(s, false)) && cmd_balanced(win_wide(s, true)),
        "Unbalanced quotes or trailing ^ in '" + s + "'");
      // Nor the command line or current directory, which cmd.exe would expand unquoted
      std::string u = s, x;
      if (win_narrow(win_wide(s, true), x)) u += "\n" + x;
      for (char& c : u) if (c>='a' && c<='z') c = static_cast<char>(c-'a'+'A');
      bool unsupported = false, cd = false;
      for (size_t i=u.find_first_of("%!"); i!=std::string::npos; i=u.find_first_of("%!", i+1)) {
        if (i+3<u.size() && u[i+1]=='=' && u[i+2]>='A' && u[i+2]<='Z' && u[i+3]==':'
            && u.find(u[i], i+4)!=std::string::npos) unsupported = true;
      }
      for (const char* v : {"%CMDCMDLINE%", "%CMDCMDLINE:", "!CMDCMDLINE!", "!CMDCMDLINE:"}) {
        if (u.find(v)!=std::string::npos) unsupported = true;
      }
      casadi_assert(!unsupported, "Unsupported variable in '" + s + "'");
      for (const char* v : {"%CD%", "%CD:", "!CD!", "!CD:",
          "%__CD__%", "%__CD__:", "!__CD__!", "!__CD__:"}) {
        if (u.find(v)!=std::string::npos) cd = true;
      }
      if (cd) {
        std::vector<wchar_t> buf(GetCurrentDirectoryW(0, nullptr) + 1);
        DWORD n = GetCurrentDirectoryW(static_cast<DWORD>(buf.size()), buf.data());
        std::wstring cwd(buf.data(), n<buf.size() ? n : 0);
        // Like CasADi's paths, unchanged for the program: no best fit, no trail bytes
        std::string a;
        bool ascii = true, ok = !cwd.empty();
        for (wchar_t c : cwd) {
          if (c>=0x80) ascii = false;
          if (c<0x20) ok = false;
        }
        casadi_assert(ok && cwd.find_first_of(L"&|<>^()%!\"")==std::wstring::npos
          && (ascii || (win_narrow(cwd, a) && win_bytes_safe())),
          "Unsupported character in the current directory for '" + s + "'");
      }
    }
#endif // _WIN32

    // Construct the compiler command
    std::stringstream cccmd;
    cccmd << compiler;
    for (auto i=compiler_flags.begin(); i!=compiler_flags.end(); ++i) {
      cccmd << " " << *i;
    }
    // Include directories (OS-agnostic: prefix with the compiler's include flag)
    for (const std::string& dir : include_dirs) {
      cccmd << " " << compiler_include_flag << quote(include_dir(dir));
    }
    cccmd << " " << compiler_setup;

    // C/C++ source file
    cccmd << " " << quote(name_);

    // Temporary object file
    cccmd << " " + compiler_output_flag << quote(obj_name_);

    // Compile into an object
    if (verbose_) casadi_message("calling \"" + cccmd.str() + "\"");
    if (run(cccmd.str())) {
      casadi_error("Compilation failed. Tried \"" + cccmd.str() + "\"");
    }

    // Link step
    std::stringstream ldcmd;
    ldcmd << linker;

    // Temporary file
    ldcmd << " " << quote(obj_name_) << " " + linker_output_flag + quote(bin_name_);

    // Add flags
    for (auto i=linker_flags.begin(); i!=linker_flags.end(); ++i) {
      ldcmd << " " << *i;
    }
    ldcmd << " " << linker_setup;

    // Compile into a shared library
    if (verbose_) casadi_message("calling \"" + ldcmd.str() + "\"");
    if (run(ldcmd.str())) {
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
