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
#include <algorithm>
#include <sstream>
#else // _WIN32
#include <cerrno>
#include <spawn.h>
#include <sys/wait.h>
#ifdef __APPLE__
#include <crt_externs.h>
#else // __APPLE__
extern char **environ;
#endif // __APPLE__
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

    // Command interpreter, for batch files
    std::string win_cmd() {
      char dir[MAX_PATH];
      UINT n = GetSystemDirectoryA(dir, MAX_PATH);
      return std::string(dir, n<MAX_PATH ? n : 0) + "\\cmd.exe";
    }

#endif // _WIN32

    // Arguments of an option, split at whitespace; quotes group, \" is a quote
    void add_args(std::vector<std::string>& args, const std::string& s) {
#ifdef _WIN32
      const char* quotes = "\"";
#else // _WIN32
      const char* quotes = "\"'";
#endif // _WIN32
      std::string a;
      bool arg = false;
      char q = 0;
      for (size_t i=0; i<s.size(); ++i) {
        char c = s[i];
        if (c=='\\' && q!='\'' && i+1<s.size() && s[i+1]=='"') {
          a += s[++i];
          arg = true;
        } else if (q ? c==q : c && std::strchr(quotes, c)) {
          q = q ? 0 : c;
          arg = true;
        } else if (!q && (c==' ' || c=='\t' || c=='\n' || c=='\r')) {
          if (arg) args.push_back(a);
          a.clear();
          arg = false;
        } else {
          a += c;
          arg = true;
        }
      }
      casadi_assert(!q, "Unbalanced quotes in '" + s + "'");
      if (arg) args.push_back(a);
    }

    // A path as one argument, not to be taken for an option
    std::string path_arg(const std::string& s) {
      if (!s.empty() && (s[0]=='-' || s[0]=='@')) return "./" + s;
      return s;
    }

    // A flag followed by a path, in one argument unless the flag ends with whitespace
    void add_path(std::vector<std::string>& args, const std::string& flag,
        const std::string& path) {
      size_t n = args.size();
      add_args(args, flag);
      if (args.size()>n && flag.find_last_of(" \t")!=flag.size()-1) {
        args.back() += path_arg(path);
      } else {
        args.push_back(path_arg(path));
      }
    }

    // Command line of arguments, as split by the MSVC runtime
    std::string cmdline(const std::vector<std::string>& args) {
      std::string r;
      for (const std::string& a : args) {
        if (!r.empty()) r += ' ';
        if (!a.empty() && a.find_first_of(" \t\"&|<>()")==std::string::npos) {
          r += a;
          continue;
        }
        r += '"';
        size_t n = 0;
        for (char c : a) {
          if (c=='\\') {
            n++;
          } else {
            r += std::string(c=='"' ? 2*n+1 : n, '\\') + c;
            n = 0;
          }
        }
        r += std::string(2*n, '\\') + '"';
      }
      return r;
    }

#ifdef _WIN32
    // Extensions cmd.exe tries
    std::vector<std::string> win_pathext() {
      const char* pathext = getenv("PATHEXT");
      std::stringstream ss(pathext ? pathext : ".COM;.EXE;.BAT;.CMD");
      std::vector<std::string> r;
      std::string e;
      while (std::getline(ss, e, ';')) if (!e.empty()) r.push_back(e);
      return r;
    }

    // Full path of a program, found as by cmd.exe on the PATH; empty if none
    std::string win_find(const std::string& prog) {
      std::vector<wchar_t> buf(32768);
      DWORD n = GetEnvironmentVariableW(L"PATH", buf.data(), static_cast<DWORD>(buf.size()));
      std::wstring path(buf.data(), n<buf.size() ? n : 0);
      std::vector<std::string> names;
      for (const std::string& e : win_pathext()) {
        if (prog.size()>e.size() && _stricmp(prog.c_str()+prog.size()-e.size(), e.c_str())==0) {
          names.insert(names.begin(), prog);
        }
        names.push_back(prog + e);
      }
      bool utf8 = GetACP()==CP_UTF8;
      for (const std::string& name : names) {
        wchar_t wname[MAX_PATH], wfull[MAX_PATH];
        if (!MultiByteToWideChar(CP_ACP, 0, name.c_str(), -1, wname, MAX_PATH)) continue;
        DWORD n = SearchPathW(path.c_str(), wname, nullptr, MAX_PATH, wfull, nullptr);
        if (n==0 || n>=MAX_PATH) continue;
        char full[4*MAX_PATH];
        BOOL lossy = FALSE;
        if (WideCharToMultiByte(CP_ACP, utf8 ? 0 : WC_NO_BEST_FIT_CHARS, wfull, -1, full,
            sizeof(full), nullptr, utf8 ? nullptr : &lossy)>0 && !lossy) return full;
      }
      return "";
    }

    // Run a program
    int run(const std::vector<std::string>& args) {
      std::vector<std::string> a = args;
      if (a.empty() || (a[0] = win_find(a[0])).empty()) return -1;
      // A batch file is run by cmd.exe, which must find the arguments within quotes
      std::string ext = a[0].substr(a[0].size()-4);
      if (_stricmp(ext.c_str(), ".bat") && _stricmp(ext.c_str(), ".cmd")) {
        return win_run(a[0], cmdline(a));
      }
      for (const std::string& s : a) {
        bool quote = s.find('"')!=std::string::npos;
        bool ok = s.find_first_of(quote ? "%!^&|<>()" : "%!^")==std::string::npos
          && std::count(s.begin(), s.end(), '"')%2==0;
        for (char c : s) ok = ok && static_cast<unsigned char>(c)>=0x20;
        casadi_assert(ok, "Unsupported argument '" + s + "' for a batch file");
      }
      std::string shell = win_cmd();
      return win_run(shell, "\"" + shell + "\" /d /s /c \"" + cmdline(a) + "\"");
    }
#else // _WIN32
    // Run a program
    int run(const std::vector<std::string>& args) {
      if (args.empty()) return -1;
      std::vector<char*> argv;
      for (const std::string& a : args) argv.push_back(const_cast<char*>(a.c_str()));
      argv.push_back(nullptr);
#ifdef __APPLE__
      char** env = *_NSGetEnviron();
#else // __APPLE__
      char** env = environ;
#endif // __APPLE__
      pid_t pid;
      int status;
      if (posix_spawnp(&pid, argv[0], nullptr, nullptr, argv.data(), env)) return -1;
      while (waitpid(pid, &status, 0)<0) if (errno!=EINTR) return -1;
      return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    }
#endif // _WIN32
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

    // Construct the compiler command
    std::vector<std::string> cccmd;
    add_args(cccmd, compiler);
    for (const std::string& f : compiler_flags) add_args(cccmd, f);
    // Include directories (OS-agnostic: prefix with the compiler's include flag)
    for (const std::string& dir : include_dirs) add_path(cccmd, compiler_include_flag, dir);
    add_args(cccmd, compiler_setup);

    // C/C++ source file
    cccmd.push_back(path_arg(name_));

    // Temporary object file
    add_path(cccmd, compiler_output_flag, obj_name_);

    // Compile into an object
    if (verbose_) casadi_message("calling \"" + cmdline(cccmd) + "\"");
    int ret = run(cccmd);
    if (ret) {
      casadi_error("Compilation failed. Tried \"" + cmdline(cccmd) + "\"");
    }

    // Link step
    std::vector<std::string> ldcmd;
    add_args(ldcmd, linker);

    // Temporary file
    ldcmd.push_back(path_arg(obj_name_));
    add_path(ldcmd, linker_output_flag, bin_name_);

    // Add flags
    for (const std::string& f : linker_flags) add_args(ldcmd, f);
    add_args(ldcmd, linker_setup);

    // Compile into a shared library
    if (verbose_) casadi_message("calling \"" + cmdline(ldcmd) + "\"");
    ret = run(ldcmd);
    if (ret) {
      casadi_error("Linking failed. Tried \"" + cmdline(ldcmd) + "\"");
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
