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


#include "modelica_parser_internal.hpp"

namespace casadi {

ModelicaParserInternal::ModelicaParserInternal() {
}

ModelicaParserInternal::~ModelicaParserInternal() {
}

std::map<std::string, ModelicaParserInternal::Plugin> ModelicaParserInternal::solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::mutex ModelicaParserInternal::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

const std::string ModelicaParserInternal::infix_ = "modelica";

void ModelicaParserInternal::disp(std::ostream &stream, bool more) const {
}

void ModelicaParserInternal::parse(const std::string& filename) {
  casadi_error("parse not defined for " + class_name());
}

} // namespace casadi
