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


#ifndef CASADI_EXCEPTION_HPP
#define CASADI_EXCEPTION_HPP

#include <exception>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include <casadi/core/casadi_export.h>

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
*/
class CasadiException : public std::exception {
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

  //! \brief Append a message
  CasadiException& operator<<(const std::string& msg) {
    msg_ += msg;
    return *this;
  }

  //! \brief Append an exception
  CasadiException& operator<<(const std::exception& ex) {
    msg_ += " => ";
    msg_ += ex.what();
    return *this;
  }

  protected:
  std::string msg_;
};

class KeyboardInterruptException : public CasadiException {
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

// Convert to string
#define CASADI_STR1(x) #x
#define CASADI_STR(x) CASADI_STR1(x)

// String denoting where the assertion is situated
#define CASADI_WHERE casadi::trim_path(__FILE__ ":" CASADI_STR(__LINE__))

// Throw an exception with information about source code location
#define casadi_error(msg) \
  throw casadi::CasadiException(CASADI_WHERE + ": " + std::string(msg));

// This assertion checks for illegal user inputs
#define casadi_assert_message(x, msg) \
  if (!(x)) casadi_error("Assertion \"" CASADI_STR(x) "\" failed:\n" + std::string(msg));

// This assertion if for errors caused by bugs in CasADi, use it instead of C:s assert(),
// but never in destructors
#if NDEBUG
#define casadi_assert(x) casadi_assert_message(x, "Please notify the CasADi developers.")
#else
#define casadi_assert(x) casadi_assert_message(x, "Please notify the CasADi developers.\n"\
    "Hint for developers: GlobalOptions.setCatchErrorsSwig(False)" \
    " to obtain gdb stacktrace in python.")
#endif

// Issue a warning, including location in the source code
#define casadi_warning(msg) \
  casadi::userOut<true, casadi::PL_WARN>() \
    << "CasADi warning at " << CASADI_WHERE << ": " << msg << "\n";

// Issue a warning if an assertion fails
#define casadi_assert_warning(x, msg) \
  if (!(x)) casadi_warning("Assertion \"" CASADI_STR(x) "\" failed.");

// Issue a message, including location in the source code
#define casadi_message(msg) \
  casadi::userOut() \
    << "CasADi message at " << CASADI_WHERE << ": " << msg << "\n";

} // namespace casadi

#endif // CASADI_EXCEPTION_HPP
