/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "casadi_os.hpp"
#include "exception.hpp"
#include "global_options.hpp"
#include <bitset>

#ifndef _WIN32
#ifdef WITH_DEEPBIND
#ifndef __APPLE__
#if __GLIBC__
extern char **environ;
#endif
#endif
#endif
#endif


#ifdef _WIN32
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#endif

namespace casadi {

// http://stackoverflow.com/questions/303562/c-format-macro-inline-ostringstream
#define STRING(ITEMS) \
  ((dynamic_cast<std::ostringstream &>(std::ostringstream() \
    . seekp(0, std::ios_base::cur) << ITEMS)) . str())

char pathsep() {
    #ifdef _WIN32
        return ';';
    #else
        return ':';
    #endif
}
std::string filesep() {
    #ifdef _WIN32
        return "\\";
    #else
        return "/";
    #endif
}

std::vector<std::string> get_search_paths() {

    // Build up search paths;
    std::vector<std::string> search_paths;

    // Search path: global casadipath option
    std::stringstream casadipaths(GlobalOptions::getCasadiPath());
    std::string casadipath;
    while (std::getline(casadipaths, casadipath, pathsep())) {
        search_paths.push_back(casadipath);
    }

    // Search path: CASADIPATH env variable
    char* pLIBDIR;
    pLIBDIR = getenv("CASADIPATH");

    if (pLIBDIR!=nullptr) {
        std::stringstream casadipaths(pLIBDIR);
        std::string casadipath;
        while (std::getline(casadipaths, casadipath, pathsep())) {
        search_paths.push_back(casadipath);
        }
    }

    // Search path: bare
    search_paths.push_back("");

    // Search path : PLUGIN_EXTRA_SEARCH_PATH
    #ifdef PLUGIN_EXTRA_SEARCH_PATH
    search_paths.push_back(
        std::string("") + PLUGIN_EXTRA_SEARCH_PATH);
    #endif // PLUGIN_EXTRA_SEARCH_PATH

    // Search path : current directory
    search_paths.push_back(".");

    return search_paths;
}

#ifdef WITH_DL

handle_t open_shared_library(const std::string& lib, const std::vector<std::string> &search_paths,
    const std::string& caller, bool global) {
        std::string resultpath;
        return open_shared_library(lib, search_paths, resultpath, caller, global);
}

int close_shared_library(handle_t handle) {
    #ifdef _WIN32
        return !FreeLibrary(handle);
    #else // _WIN32
        return dlclose(handle);
    #endif // _WIN32
}

handle_t open_shared_library(const std::string& lib, const std::vector<std::string> &search_paths,
        std::string& resultpath, const std::string& caller, bool global) {
    // Alocate a handle
    handle_t handle = 0;

    // Alocate a handle pointer
    #ifndef _WIN32
        int flag;
        if (global) {
            flag = RTLD_NOW | RTLD_GLOBAL;
        } else {
            flag = RTLD_LAZY | RTLD_LOCAL;
        }
    #ifdef WITH_DEEPBIND
    #if !defined(__APPLE__) && !defined(__EMSCRIPTEN__)
        flag |= RTLD_DEEPBIND;

        #if __GLIBC__
        // Workaround for https://github.com/conda-forge/casadi-feedstock/issues/93
        // and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=111556
        // In a nutshell, if RTLD_DEEPBIND is used and multiple symbols of environ
        // (one in executable's .bss and one in glibc .bss)
        // are present in the process due to copy relocations, make sure that the
        // environ in glibc .bss has the same value of environ in executable .bss
        // To avoid that over time the two values diverse due to the use of setenv,
        // we restore the original value of glibc .bss's environ at the end of the function

        // Check if there is a duplicate environ
        char*** p_environ_rtdl_next = reinterpret_cast<char ***>(dlsym(RTLD_NEXT, "environ"));
        bool environ_rtdl_next_overridden = false;
        char** environ_rtld_next_original_value = NULL;
        if (p_environ_rtdl_next && p_environ_rtdl_next != &environ) {
          environ_rtld_next_original_value = *p_environ_rtdl_next;
          *p_environ_rtdl_next = environ;
          environ_rtdl_next_overridden = true;
        }
        #endif
    #endif
    #endif
    #endif


    // Prepare error string
    std::stringstream errors;
    errors << caller << ": Cannot load shared library '"
           << lib << "': " << std::endl;
    errors << "   (\n"
           << "    Searched directories: 1. casadipath from GlobalOptions\n"
           << "                          2. CASADIPATH env var\n"
           << "                          3. PATH env var (Windows)\n"
           << "                          4. LD_LIBRARY_PATH env var (Linux)\n"
           << "                          5. DYLD_LIBRARY_PATH env var (osx)\n"
           << "    A library may be 'not found' even if the file exists:\n"
           << "          * library is not compatible (different compiler/bitness)\n"
           << "          * the dependencies are not found\n"
           << "   )";

    std::string searchpath;

    // Try getting a handle
    for (casadi_int i=0;i<search_paths.size();++i) {
      searchpath = search_paths[i];
#ifdef _WIN32
      SetDllDirectory(TEXT(searchpath.c_str()));
      handle = LoadLibrary(TEXT(lib.c_str()));
      SetDllDirectory(NULL);
#else // _WIN32
      std::string libname = searchpath.empty() ? lib : searchpath + filesep() + lib;
      handle = dlopen(libname.c_str(), flag);
#endif // _WIN32
      if (handle) {
        resultpath = searchpath;
        break;
      } else {
        errors << std::endl << "  Tried '" << searchpath << "' :";
#ifdef _WIN32
        errors << std::endl << "    Error code (WIN32): " << STRING(GetLastError());
#else // _WIN32
        errors << std::endl << "    Error code: " << dlerror();
#endif // _WIN32
      }
    }

    #ifndef _WIN32
    #ifdef WITH_DEEPBIND
    #ifndef __APPLE__
    #if __GLIBC__
        if (environ_rtdl_next_overridden) {
          *p_environ_rtdl_next = environ_rtld_next_original_value;
          environ_rtdl_next_overridden = false;
        }
    #endif
    #endif
    #endif
    #endif

    casadi_assert(handle!=nullptr, errors.str());

    return handle;
}

#endif // WITH_DL


// Convert UTF-8 to UTF-16 on Windows
#ifdef _WIN32
std::wstring utf8_to_utf16(const std::string& s) {
    int wlen = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
    if (wlen == 0) return {};
    std::wstring ws(wlen, 0);
    MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), &ws[0], wlen);
    return ws;
}
#endif

#ifdef _WIN32
class FdStreamBuf : public std::streambuf {
public:
    explicit FdStreamBuf(int fd, size_t bufsize = 4096)
        : fd_(fd), buffer_(bufsize) {
        setg(buffer_.data(), buffer_.data(), buffer_.data());
    }

    ~FdStreamBuf() override {
        if (fd_ >= 0) {
            _close(fd_);
        }
    }

protected:
    int_type underflow() override {
        if (gptr() < egptr()) {
            return traits_type::to_int_type(*gptr());
        }

        int n = _read(fd_, buffer_.data(), static_cast<unsigned int>(buffer_.size()));
        if (n <= 0) {
            return traits_type::eof();
        }

        setg(buffer_.data(), buffer_.data(), buffer_.data() + n);
        return traits_type::to_int_type(*gptr());
    }

    std::streampos seekoff(std::streamoff off, std::ios_base::seekdir dir,
                            std::ios_base::openmode which = std::ios_base::in) override {
        if (!(which & std::ios_base::in)) return -1;

        __int64 whence;
        switch (dir) {
            case std::ios_base::beg: whence = SEEK_SET; break;
            case std::ios_base::cur:
                // Need to include the offset in buffer
                off -= egptr() - gptr();
                whence = SEEK_CUR;
                break;
            case std::ios_base::end: whence = SEEK_END; break;
            default: return -1;
        }

        __int64 result = _lseeki64(fd_, off, static_cast<int>(whence));
        if (result == -1) {
            return -1;
        }

        // Invalidate the buffer
        setg(buffer_.data(), buffer_.data(), buffer_.data());
        return result;
    }

    std::streampos seekpos(std::streampos pos,
                            std::ios_base::openmode which = std::ios_base::in) override {
        return seekoff(static_cast<std::streamoff>(pos), std::ios_base::beg, which);
    }


private:
    int fd_;
    std::vector<char> buffer_;
};

class FdOStreamBuf : public std::streambuf {
public:
    explicit FdOStreamBuf(int fd, size_t bufsize = 4096)
        : fd_(fd), buffer_(bufsize) {
        setp(buffer_.data(), buffer_.data() + buffer_.size());
    }

    ~FdOStreamBuf() override {
        sync(); // flush on destruction
        if (fd_ >= 0) {
            _close(fd_);
        }
    }

protected:
    int_type overflow(int_type ch) override {
        if (flush_buffer() == -1) return traits_type::eof();

        if (ch != traits_type::eof()) {
            *pptr() = static_cast<char>(ch);
            pbump(1);
        }

        return ch;
    }

    int sync() override {
        return flush_buffer() == -1 ? -1 : 0;
    }

private:
    int flush_buffer() {
        int len = static_cast<int>(pptr() - pbase());
        if (len > 0) {
            int written = _write(fd_, pbase(), len);
            if (written != len) return -1;
            pbump(-len);
        }
        return 0;
    }

    int fd_;
    std::vector<char> buffer_;
};

struct OwnedIStream {
    std::unique_ptr<FdStreamBuf> buffer;
    std::unique_ptr<std::istream> stream;

    explicit OwnedIStream(int fd)
        : buffer(std::make_unique<FdStreamBuf>(fd)),
          stream(std::make_unique<std::istream>(buffer.get())) {}
};

struct StreamWithOwnedBuffer : public std::istream {
    std::shared_ptr<OwnedIStream> owned;
    explicit StreamWithOwnedBuffer(std::shared_ptr<OwnedIStream> o)
        : std::istream(o->buffer.get()), owned(std::move(o)) {}
};

struct OwnedOStream {
    std::unique_ptr<FdOStreamBuf> buffer;
    std::unique_ptr<std::ostream> stream;

    explicit OwnedOStream(int fd)
        : buffer(std::make_unique<FdOStreamBuf>(fd)),
          stream(std::make_unique<std::ostream>(buffer.get())) {}
};

struct StreamWithOwnedOBuffer : public std::ostream {
    std::shared_ptr<OwnedOStream> owned;
    explicit StreamWithOwnedOBuffer(std::shared_ptr<OwnedOStream> o)
        : std::ostream(o->buffer.get()), owned(std::move(o)) {}
};
#endif

// Portable ifstream opener that supports UTF-8 filenames on Windows
std::unique_ptr<std::istream> ifstream_compat(const std::string& utf8_path,
                                    std::ios::openmode mode) {
#ifdef _WIN32
    std::wstring utf16_path = utf8_to_utf16(utf8_path);
    DWORD access = 0;
    access |= GENERIC_READ;
    if (mode & std::ios::out) access |= GENERIC_WRITE;

    HANDLE h = CreateFileW(utf16_path.c_str(), access,
                             FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr,
                           OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (h == INVALID_HANDLE_VALUE) return {};

    int flags = (mode & std::ios::out) ? _O_RDWR : _O_RDONLY;
    if (mode & std::ios::binary) flags |= _O_BINARY;

    int fd = _open_osfhandle(reinterpret_cast<intptr_t>(h), flags);
    if (fd == -1) {
        CloseHandle(h);
        return {};
    }
    auto owned = std::make_shared<OwnedIStream>(fd);
    return std::unique_ptr<StreamWithOwnedBuffer>(new StreamWithOwnedBuffer(std::move(owned)));
#else
    auto ifs = std::unique_ptr<std::ifstream>(new std::ifstream(utf8_path, mode));
    if (!*ifs) return {};
    return std::unique_ptr<std::istream>(std::move(ifs));
#endif
}

std::unique_ptr<std::ostream> ofstream_compat(const std::string& utf8_path,
                                              std::ios::openmode mode) {
#ifdef _WIN32
    std::wstring utf16_path = utf8_to_utf16(utf8_path);

    DWORD access = 0;
    access |= GENERIC_WRITE;
    if (mode & std::ios::in)   access |= GENERIC_READ;

    DWORD creation = (mode & std::ios::app) ? OPEN_ALWAYS : CREATE_ALWAYS;

    HANDLE h = CreateFileW(utf16_path.c_str(), access,
                            FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr,
                           creation, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (h == INVALID_HANDLE_VALUE) return {};

    int flags = (mode & std::ios::in) ? _O_RDWR : _O_WRONLY;
    if (mode & std::ios::app) flags |= _O_APPEND;
    if (mode & std::ios::binary) flags |= _O_BINARY;

    int fd = _open_osfhandle(reinterpret_cast<intptr_t>(h), flags);
    if (fd == -1) {
        CloseHandle(h);
        return {};
    }

    auto owned = std::make_shared<OwnedOStream>(fd);
    return std::unique_ptr<StreamWithOwnedOBuffer>(new StreamWithOwnedOBuffer(std::move(owned)));
#else
    auto ofs = std::unique_ptr<std::ofstream>(new std::ofstream(utf8_path, mode));
    if (!*ofs) return {};
    return std::unique_ptr<std::ostream>(std::move(ofs));
#endif
}

} // namespace casadi
