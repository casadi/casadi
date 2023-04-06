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


#ifndef CASADI_EXCEPTION_HPP
#define CASADI_EXCEPTION_HPP

#include <chrono>
#include <ctime>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <casadi/core/casadi_export.h>

// Disable some Visual studio warnings
#ifdef _MSC_VER

// warning C4251: Need a dll interface?
#pragma warning(disable:4251)

// warning C4275: non dll-interface class 'std::exception' used as base for dll-interface
// class 'casadi::CasadiException'
#pragma warning(disable:4275)

// warning C4996: 'sprintf': This function or variable may be unsafe. Consider using sprintf_s
// instead
#pragma warning(disable:4996)

#endif // _MSC_VER

namespace casadi {

/** \brief  Casadi exception class

        \author Joel Andersson
        \date 2010
        Example for simple exception throwing:
        \code
                throw CasadiException("This is a nasty error");
        \endcode
        Example for exception chaining:
        \code
                try {
                        throw CasadiException("This is a nasty error");
                catch(CasadiException &e) {
                        throw CasadiException("Serious error.") << e;
                }
        \endcode

    \identifier{7u} */
class CASADI_EXPORT CasadiException : public std::exception {
  public:
  //! \brief Default constructor
  CasadiException() {
  }

  //! \brief Form message string
  explicit CasadiException(const std::string& msg) : msg_(msg) {}

  //! \brief Destructor
  ~CasadiException() throw() {}

  //! \brief Display error
  const char* what() const throw() override {
    return msg_.c_str();
  }

  protected:
  std::string msg_;
};

class CASADI_EXPORT KeyboardInterruptException : public CasadiException {
  public:
  //! \brief Default constructor
  KeyboardInterruptException() : CasadiException("KeyboardInterrupt") {}
  //! \brief Destructor
  ~KeyboardInterruptException() throw() {}
};

// Strip path prefix
inline std::string trim_path(const std::string& full_path) {
  size_t found = full_path.rfind("/casadi/");
  if (found == std::string::npos) {
    return full_path;
  } else {
    std::string ret = full_path;
    ret.replace(0, found, "...");
    return ret;
  }
}

// Current time as a string
inline std::ostream& message_prefix(std::ostream &stream) {
  // CasADi prefix
  stream << "CasADi - ";
  // Get current time
  auto now = std::chrono::system_clock::now();
  std::time_t tt = std::chrono::system_clock::to_time_t(now);
  auto local_tm = *std::localtime(&tt);  // NOLINT(runtime/threadsafe_fn)
  // Convert to YYYY-MM-DD HH:MM:SS format
  stream << local_tm.tm_year + 1900 << '-';  // YYYY-
  stream << std::setfill('0') << std::setw(2) << local_tm.tm_mon + 1 << '-';  // MM-
  stream << std::setfill('0') << std::setw(2) << local_tm.tm_mday << ' ';  // DD
  stream << std::setfill('0') << std::setw(2) << local_tm.tm_hour << ':';  // hh:
  stream << std::setfill('0') << std::setw(2) << local_tm.tm_min << ':';  // mm:
  stream << std::setfill('0') << std::setw(2) << local_tm.tm_sec;  // ss
  return stream;
}

// String denoting where the macro is situated
#define CASADI_WHERE casadi::trim_path(__FILE__ ":" CASADI_STR(__LINE__))

// Throw an exception with information about source code location
#define casadi_error(msg, ...) \
throw casadi::CasadiException(CASADI_WHERE + ": "\
          + casadi::fmtstr(msg, casadi::strvec(__VA_ARGS__)))

// This assertion checks for illegal user inputs
#define casadi_assert(x, msg, ...) \
if (!(x)) casadi_error("Assertion \"" CASADI_STR(x) "\" failed:\n"\
          + std::string(msg), __VA_ARGS__)

// This assertion if for internal errors caused by bugs in CasADi
#define casadi_assert_dev(x) casadi_assert(x, "Notify the CasADi developers.")

// This assertion if for internal errors caused by bugs in CasADi
#define casadi_report() casadi_error("Notify the CasADi developers.")

// Issue a warning, including location in the source code
#define casadi_warning(msg) \
  casadi::message_prefix(casadi::uerr()) \
    << " WARNING(\"" << msg << "\") [" << CASADI_WHERE << "]\n" << std::flush;

// Issue a message, including location in the source code
#define casadi_message(msg) \
  casadi::message_prefix(casadi::uout()) \
    << " MESSAGE(\"" << msg << "\") [" << CASADI_WHERE << "]\n" << std::flush;

} // namespace casadi

#endif // CASADI_EXCEPTION_HPP
