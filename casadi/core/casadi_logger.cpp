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

#include "casadi_common.hpp"
#include "casadi_logger.hpp"

#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

namespace casadi {

#ifdef CASADI_WITH_THREAD
  std::mutex mutex_logger;
#endif //CASADI_WITH_THREAD

  void Logger::WriteFunThreadSafe(const char* s, std::streamsize num, bool error) {
#ifdef CASADI_WITH_THREAD
    std::lock_guard<std::mutex> lock(mutex_logger);
#endif //CASADI_WITH_THREAD
    writeFun(s, num, error);
  }

  void Logger::FlushThreadSafe(bool error) {
#ifdef CASADI_WITH_THREAD
    std::lock_guard<std::mutex> lock(mutex_logger);
#endif //CASADI_WITH_THREAD
    flush(error);
  }

  void (*Logger::writeFun)(const char* s, std::streamsize num, bool error) =
    Logger::writeDefault;

  void (*Logger::flush)(bool error) =Logger::flushDefault;

  std::ostream& uout() {
    // Singleton pattern
    static Logger::Stream<false> instance;
    return instance;
  }

  std::ostream& uerr() {
    // Singleton pattern
    static Logger::Stream<true> instance;
    return instance;
  }

} // namespace casadi

extern "C" {
  int casadi_printf(const char *fmt, ...) {
    // Variable number of arguments
    va_list args;
    va_start(args, fmt);
    // Static & dynamic buffers
    char buf[256];
    size_t buf_sz = sizeof(buf);
    char* buf_dyn = nullptr;
    // Try to print with a small buffer
    casadi_int n = vsnprintf(buf, buf_sz, fmt, args);
    // Need a larger buffer?
    if (n>static_cast<casadi_int>(buf_sz)) {
      buf_sz = static_cast<size_t>(n+1);
      buf_dyn = new char[buf_sz];
      n = vsnprintf(buf_dyn, buf_sz, fmt, args);
    }
    // Print buffer content
    if (n>=0) casadi::uout() << (buf_dyn ? buf_dyn : buf) << std::flush;
    // Cleanup
    delete[] buf_dyn;
    va_end(args);
    return n;
  }
}
