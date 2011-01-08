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

#include "xml_arg.hpp"
#include <casadi/sx/sx_tools.hpp>

namespace CasADi{
namespace Modelica{

StrArg::StrArg(const std::string& str) : str(str){
}

StrArg::operator std::string() const{
  return str;
}

StrArg::operator int() const{
  std::istringstream buffer(str);
  int ret;
  buffer >> ret;
  return ret;
}

StrArg::operator double() const{
  std::istringstream buffer(str);
  double ret;
  buffer >> ret;
  return ret;
}

StrArg::operator SXMatrix() const{
  std::istringstream buffer(str);
  SXMatrix ret;
  buffer >> ret;
  return ret;
}

StrArg::operator SX() const{
  SXMatrix m = *this;
  return SX(m.getElement());
}

StrArg::operator MX() const{
  return MX(str);
}


std::ostream& operator<<(std::ostream &stream, const StrArg& arg){
  return stream << arg.str;
}

} // namespace Modelica
} // namespace CasADi

