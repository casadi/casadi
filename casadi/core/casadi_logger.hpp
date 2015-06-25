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


#ifndef CASADI_LOGGER_HPP
#define CASADI_LOGGER_HPP

#include "casadi_common.hpp"

#include <iostream>
#include <fstream>

namespace casadi {
  /**
   * \brief Keeps track of logging output to screen and/or files.
   * All printout from CasADi routines should go through this files.
   *
   *  \author Joel Andersson
   *  \date 2015
   */
  class CASADI_EXPORT Logger {
  private:
    /// No implementation - no instances are allowed of this class
    Logger();
  public:
    /// Print output message, can be redefined
    static void (*writeOut)(const char* s, std::streamsize num);

    /// Print error message, can be redefined
    static void (*writeErr)(const char* s, std::streamsize num);

    /// Print output message, default
    static void writeOutDefault(const char* s, std::streamsize num) {
      std::cout.write(s, num);
    }

    /// Print error message, default
    static void writeErrDefault(const char* s, std::streamsize num) {
      std::cerr.write(s, num);
    }

    /// Print log message, single character
    static void writeOutCh(char ch) { writeOut(&ch, 1);}

    /// Print error message, single character
    static void writeErrCh(char ch) { writeOut(&ch, 1);}
  };

  // Stream buffer for csout like printing
  template<bool Err>
  class LoggerStreambuf : public std::streambuf {
  public:
    LoggerStreambuf() {}
  protected:
    virtual int_type overflow(int_type ch) {
      if (ch != traits_type::eof()) {
        if (Err) {
          Logger::writeErrCh(static_cast<char>(ch));
        } else {
          Logger::writeOutCh(static_cast<char>(ch));
        }
      }
      return ch;
    }
    virtual std::streamsize xsputn(const char* s, std::streamsize num) {
      if (Err) {
        Logger::writeErr(s, num);
      } else {
        Logger::writeOut(s, num);
      }
      return num;
    }
  };

  // Output stream for csout like printing
  template<bool Err>
  class LoggerStream : public std::ostream {
    protected:
    LoggerStreambuf<Err> buf;
  public:
    LoggerStream() : std::ostream(&buf) {}
  };

  // This replaces csout
  static LoggerStream<false> csout;

  // This replaces cserr
  static LoggerStream<true> cserr;

} // namespace casadi

#endif // CASADI_LOGGER_HPP
