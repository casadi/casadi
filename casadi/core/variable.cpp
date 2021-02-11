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


#include "exception.hpp"
#include "variable.hpp"

using namespace std;
namespace casadi {

const char* to_string(Causality v) {
  switch (v) {
  case PARAMETER: return "parameter";
  case CALCULATED_PARAMETER: return "calculatedParameter";
  case INPUT: return "input";
  case OUTPUT: return "output";
  case LOCAL: return "local";
  case INDEPENDENT: return "independent";
  }
  return 0;
}

const char* to_string(Variability v) {
  switch (v) {
  case CONSTANT: return "constant";
  case FIXED: return "fixed";
  case TUNABLE: return "tunable";
  case DISCRETE: return "discrete";
  case CONTINUOUS: return "continuous";
  }
  return 0;
}

Variable::Variable(const std::string& name, const Sparsity& sp, const MX& v, const MX& d)
    : v(v), d(d) {
  if (this->v.is_empty()) this->v = MX::sym(name, sp);
  if (this->d.is_empty()) this->d = MX::sym("der_" + name, sp);
  this->variability = CONTINUOUS;
  this->causality = LOCAL;
  this->category = CAT_UNKNOWN;
  this->description = "";
  this->valueReference = -1;
  this->min = -numeric_limits<double>::infinity();
  this->max = numeric_limits<double>::infinity();
  this->guess = 0;
  this->nominal = 1.0;
  this->start = 0.0;
  this->derivative_start = 0.0;
  this->unit = "";
  this->display_unit = "";
  this->free = false;
}

string Variable::name() const {
  return this->v.name();
}

void Variable::disp(ostream &stream, bool more) const {
  stream << name();
}

} // namespace casadi
