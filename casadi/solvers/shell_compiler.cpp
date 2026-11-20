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
#include <map>
#include <sstream>
#ifdef CASADI_WITH_THREAD
#include <mutex>
#endif // CASADI_WITH_THREAD
#else // _WIN32
#include <unistd.h>
#endif // _WIN32

namespace casadi {

namespace {
#ifdef _WIN32
    // Variables set up for the compiler and linker
    struct WinEnv {
      std::wstring block;
      std::wstring path;
    };

    // Case-insensitive order of environment variable names
    struct WinLess {
      bool operator()(const std::wstring& a, const std::wstring& b) const {
        return _wcsicmp(a.c_str(), b.c_str())<0;
      }
    };

    // Current environment with variables added or changed
    std::wstring win_block(const std::wstring& changes) {
      std::map<std::wstring, std::wstring, WinLess> vars;
      wchar_t* cur = GetEnvironmentStringsW();
      for (const wchar_t* p : {static_cast<const wchar_t*>(cur), changes.c_str()}) {
        for (; p && *p; p += wcslen(p)+1) {
          std::wstring v = p;
          size_t eq = v.find(L'=', 1);
          if (eq!=std::wstring::npos) vars[v.substr(0, eq)] = v.substr(eq+1);
        }
      }
      if (cur) FreeEnvironmentStringsW(cur);
      std::wstring r;
      for (auto&& kv : vars) {
        r += kv.first + L"=" + kv.second;
        r += L'\0';
      }
      return r + L'\0';
    }

    // Run a command line, optionally capturing stdout; -1 if it could not be started
    int win_run(const std::string& app, const std::string& cmd, const WinEnv* env,
        std::string* out=nullptr, const char* dir=nullptr) {
      HANDLE rd = nullptr, wr = nullptr;
      STARTUPINFOA si;
      ZeroMemory(&si, sizeof(si));
      si.cb = sizeof(si);
      if (out) {
        if (!CreatePipe(&rd, &wr, nullptr, 0)) return -1;
        si.dwFlags = STARTF_USESTDHANDLES;
        si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
        si.hStdOutput = wr;
        si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
      }
      // No console window when the parent has none
      DWORD flags = GetConsoleCP() ? 0 : CREATE_NO_WINDOW;
      if (env) flags |= CREATE_UNICODE_ENVIRONMENT;
      std::vector<char> buf(cmd.begin(), cmd.end());
      buf.push_back('\0');
      std::wstring block = env ? win_block(env->block) : L"";
      PROCESS_INFORMATION pi;
      BOOL ok;
      {
#ifdef CASADI_WITH_THREAD
        // Only this child inherits the pipe
        static std::mutex mtx;
        std::lock_guard<std::mutex> lock(mtx);
#endif // CASADI_WITH_THREAD
        if (wr) SetHandleInformation(wr, HANDLE_FLAG_INHERIT, HANDLE_FLAG_INHERIT);
        ok = CreateProcessA(app.c_str(), buf.data(), nullptr, nullptr, TRUE, flags,
          env ? const_cast<wchar_t*>(block.data()) : nullptr, dir, &si, &pi);
        if (wr) CloseHandle(wr);
      }
      if (!ok) {
        if (rd) CloseHandle(rd);
        return -1;
      }
      if (out) {
        char chunk[4096];
        DWORD n;
        while (ReadFile(rd, chunk, sizeof(chunk), &n, nullptr) && n>0) out->append(chunk, n);
        CloseHandle(rd);
      }
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
      // POSIX shell
      if (!s.empty() && s.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-+=/.,:@%")==std::string::npos) return s;
      std::string r = "'";
      for (char c : s) {
        if (c=='\'') {
          r += "'\\''";
        } else {
          r += c;
        }
      }
      return r + "'";
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
      if (s.size()>=2 && s.front()=='\'' && s.back()=='\'') return s.substr(1, s.size()-2);
      std::string e = s;
      if (e.size()>=2 && e.front()=='"' && e.back()=='"') e = e.substr(1, e.size()-2);
      // $VAR and ${VAR}
      std::string r;
      for (size_t i=0; i<e.size(); ++i) {
        if (e[i]=='$') {
          bool brace = i+1<e.size() && e[i+1]=='{';
          size_t j = brace ? i+2 : i+1, k = j;
          while (k<e.size() && (e[k]=='_' || (e[k]>='0' && e[k]<='9')
              || (e[k]>='A' && e[k]<='Z') || (e[k]>='a' && e[k]<='z'))) k++;
          if (k>j && (!brace || (k<e.size() && e[k]=='}'))) {
            std::string name = e.substr(j, k-j);
            char buf[4096];
            const char* v = name=="PWD" ? getcwd(buf, sizeof(buf)) : nullptr;
            if (!v) v = getenv(name.c_str());
            if (v) r += v;
            i = brace ? k : k-1;
            continue;
          }
        }
        r += e[i];
      }
      return r;
#endif // _WIN32
    }

#ifdef _WIN32
    // Program of a command line, found as argv[0] by the MSVC runtime
    std::string program(const std::string& cmd) {
      size_t b = cmd.find_first_not_of(" \t");
      if (b==std::string::npos) return "";
      if (cmd[b]=='"') return cmd.substr(b+1, cmd.find('"', b+1)-b-1);
      return cmd.substr(b, cmd.find_first_of(" \t", b)-b);
    }

    // Extensions cmd.exe tries
    std::vector<std::string> win_pathext() {
      const char* pathext = getenv("PATHEXT");
      std::stringstream ss(pathext ? pathext : ".COM;.EXE;.BAT;.CMD");
      std::vector<std::string> r;
      std::string e;
      while (std::getline(ss, e, ';')) if (!e.empty()) r.push_back(e);
      return r;
    }

    // Find a program on the PATH of an environment
    bool win_find(const std::string& prog, const WinEnv* env, const std::string& ext) {
      wchar_t wprog[MAX_PATH], wext[MAX_PATH], wfull[MAX_PATH];
      if (!MultiByteToWideChar(CP_ACP, 0, prog.c_str(), -1, wprog, MAX_PATH)) return false;
      if (!MultiByteToWideChar(CP_ACP, 0, ext.c_str(), -1, wext, MAX_PATH)) return false;
      DWORD n = SearchPathW(env ? env->path.c_str() : nullptr, wprog, wext,
        MAX_PATH, wfull, nullptr);
      return n>0 && n<MAX_PATH;
    }

    // Environment after running a batch file, cached; null on failure
    const WinEnv* win_env(const std::string& bat, const std::string& args) {
      static std::map<std::string, WinEnv> cache;
#ifdef CASADI_WITH_THREAD
      static std::mutex mtx;
      std::lock_guard<std::mutex> lock(mtx);
#endif // CASADI_WITH_THREAD
      std::string key = bat + " " + args;
      auto it = cache.find(key);
      if (it==cache.end()) {
        std::string qbat;
        try {
          qbat = quote(bat);
        } catch (std::exception&) {
          return nullptr;
        }
        std::string shell = win_shell();
        std::string cmd = "\"" + shell
          + "\" /u /d /s /c \"set VSCMD_SKIP_SENDTELEMETRY=1&& " + qbat + " " + args
          + " >nul && set\"";
        std::string dir = bat.substr(0, bat.find_last_of('\\')), out;
        if (win_run(shell, cmd, nullptr, &out, dir.c_str())) out.clear();
        // Output is UTF-16, one NAME=value per line
        std::wstring wout(reinterpret_cast<const wchar_t*>(out.data()),
          out.size()/sizeof(wchar_t));
        WinEnv env;
        size_t pos = 0;
        while (pos<wout.size()) {
          size_t end = wout.find(L'\n', pos);
          if (end==std::wstring::npos) end = wout.size();
          std::wstring line = wout.substr(pos, end-pos);
          pos = end+1;
          if (!line.empty() && line.back()==L'\r') line.pop_back();
          size_t eq = line.find(L'=');
          if (eq==std::wstring::npos || eq==0) continue;
          // Only what the batch file changed
          std::wstring name = line.substr(0, eq), value = line.substr(eq+1);
          std::vector<wchar_t> cur(32768);
          DWORD n = GetEnvironmentVariableW(name.c_str(), cur.data(),
            static_cast<DWORD>(cur.size()));
          if (n && n<cur.size() && value==cur.data()) continue;
          if (_wcsicmp(name.c_str(), L"PATH")==0) env.path = value;
          env.block += line;
          env.block += L'\0';
        }
        if (!env.block.empty()) env.block += L'\0';
        it = cache.insert({key, env}).first;
      }
      return it->second.block.empty() ? nullptr : &it->second;
    }

    // Locate vcvarsall.bat of the latest Visual Studio with a C++ toolset
    std::string win_vcvarsall() {
      const char* pf = getenv("ProgramFiles(x86)");
      if (!pf) pf = getenv("ProgramFiles");
      if (!pf) return "";
      std::string dir = std::string(pf) + "\\Microsoft Visual Studio\\Installer";
#if defined(_M_ARM64) || defined(__aarch64__)
      std::string component = "Microsoft.VisualStudio.Component.VC.Tools.ARM64";
#else
      std::string component = "Microsoft.VisualStudio.Component.VC.Tools.x86.x64";
#endif
      std::string out, vswhere = dir + "\\vswhere.exe";
      if (win_run(vswhere, "\"" + vswhere + "\" -latest -products * -requires " + component
          + " -property installationPath", nullptr, &out, dir.c_str())) return "";
      size_t end = out.find_first_of("\r\n");
      if (end!=std::string::npos) out.resize(end);
      if (out.empty()) return "";
      return out + "\\VC\\Auxiliary\\Build\\vcvarsall.bat";
    }
#endif // _WIN32

    // Run a command through the shell
#ifdef _WIN32
    int run(const std::string& cmd, const WinEnv* env) {
      // As system(), without a console window
      std::string shell = win_shell();
      return win_run(shell, "\"" + shell + "\" /s /c \"" + cmd + "\"", env);
#else // _WIN32
    int run(const std::string& cmd) {
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

    // Set up Visual Studio if cl.exe is not on the PATH
    const WinEnv* env = nullptr;
    std::string cc;
    std::string prog = program(win_narrow(win_wide(compiler, true), cc) ? cc : compiler);
    if (_stricmp(prog.c_str(), "cl")==0 || _stricmp(prog.c_str(), "cl.exe")==0) {
      bool found = false;
      for (const std::string& e : win_pathext()) {
        if (win_find(prog, nullptr, e)) {
          found = true;
          break;
        }
      }
      if (!found) {
        // Located once per process
        static const std::string bat = win_vcvarsall();
        if (!bat.empty()) {
#if defined(_M_ARM64) || defined(__aarch64__)
          env = win_env(bat, "arm64");
#elif defined(_WIN64)
          env = win_env(bat, "x64");
#else
          env = win_env(bat, "x86");
#endif
          // Otherwise as before
          if (env && !win_find(prog, env, ".exe")) env = nullptr;
        }
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
#ifdef _WIN32
    int ret = run(cccmd.str(), env);
#else // _WIN32
    int ret = run(cccmd.str());
#endif // _WIN32
    if (ret) {
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
#ifdef _WIN32
    ret = run(ldcmd.str(), env);
#else // _WIN32
    ret = run(ldcmd.str());
#endif // _WIN32
    if (ret) {
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
