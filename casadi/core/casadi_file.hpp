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


#ifndef CASADI_FILE_HPP
#define CASADI_FILE_HPP

#include "casadi_logger.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

namespace casadi {

  /** \brief A parsed file
      \author Joel Andersson
      \date 2016
  */
  class CASADI_EXPORT ParsedFile {
  public:
    /** \brief Default constructor (no commands) */
    ParsedFile() {}

    /** \brief Construct from a file */
    explicit ParsedFile(const std::string& fname);

    /** \brief Construct from a file */
    explicit ParsedFile(const std::vector<std::string>& lines, int offset=0);

    /** \brief Parse a file */
    void parse(const std::string& fname);

    /** \brief Parse a list of strings */
    void parse(const std::vector<std::string>& lines, int offset);

    /** \brief Print parsed file */
    void print(std::ostream &stream=casadi::userOut()) const;

    /** \brief Does an entry exist? */
    bool has(const std::string& cmd) const;

    /** \brief Get entry as a text */
    std::string to_text(const std::string& cmd) const;

    /** \brief Convert to a type */
    template<typename T>
      T to(const std::string& cmd) const {
      std::istringstream ss(to_text(cmd));
      T ret;
      ss >> ret;
      return ret;
    }

    /** \brief Get entry as a string */
    std::string to_string(const std::string& cmd) const { return to<std::string>(cmd);}

    /** \brief Get entry as an integer */
    int to_int(const std::string& cmd) const { return to<int>(cmd);}

    /** \brief Map of commands */
    std::map<std::string, std::pair<int, std::string> > commands;
  };


} // namespace casadi

#endif // CASADI_FILE_HPP
