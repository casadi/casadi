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


#include "../casadi_exception.hpp"
#include "variable.hpp"

using namespace std;
namespace casadi {

  Variable::Variable(const std::string& name) {
    this->v = SXElement::sym(name);
    this->d = SXElement::sym("der_" + name);
    this->variability = CONTINUOUS;
    this->causality = INTERNAL;
    this->category = CAT_UNKNOWN;
    this->alias = NO_ALIAS;
    this->description = "";
    this->valueReference = -1;
    this->min = -numeric_limits<double>::infinity();
    this->max = numeric_limits<double>::infinity();
    this->initialGuess = 0;
    this->nominal = 1.0;
    this->start = 0.0;
    this->derivativeStart = 0.0;
    this->unit = "";
    this->displayUnit = "";
    this->free = false;
  }

  string Variable::name() const {
    return this->v.getName();
  }

  void Variable::repr(ostream &stream, bool trailing_newline) const {
    stream << "Variable(";
    print(stream, false);
    stream << ")";
    if (trailing_newline) stream << std::endl;
  }

  void Variable::print(ostream &stream, bool trailing_newline) const {
    stream << name();
    if (trailing_newline) stream << std::endl;
  }

} // namespace casadi
