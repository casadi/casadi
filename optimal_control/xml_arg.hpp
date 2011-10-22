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

#ifndef XML_ARG_HPP
#define XML_ARG_HPP

#include <string>
#include <sstream>

namespace CasADi{
namespace OptimalControl{

class StrArg{
  public:
  explicit StrArg(const std::string& str);

/** \brief  automatic typeconversion to various types */
  operator std::string() const;
  operator int() const;
  operator double() const;

/** \brief  print */
  friend std::ostream& operator<<(std::ostream &stream, const StrArg& arg);

  std::string str;
};

} // namespace OptimalControl
} // namespace CasADi

#endif //XML_ARG_HPP
