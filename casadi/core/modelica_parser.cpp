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

ModelicaParser::ModelicaParser() {
}

ModelicaParser::ModelicaParser(const std::string& name) {
  own(ModelicaParserInternal::getPlugin(name).creator());
}

ModelicaParser::~ModelicaParser() {
}

const ModelicaParserInternal* ModelicaParser::operator->() const {
  return static_cast<const ModelicaParserInternal*>(SharedObject::operator->());
}

ModelicaParserInternal* ModelicaParser::operator->() {
  return static_cast<ModelicaParserInternal*>(SharedObject::operator->());
}

void ModelicaParser::parse(const std::string& filename) {
  return (*this)->parse(filename);
}

void ModelicaParser::load_plugin(const std::string& name) {
  ModelicaParserInternal::load_plugin(name);
}

std::string ModelicaParser::doc(const std::string& name) {
  return ModelicaParserInternal::getPlugin(name).doc;
}

} // namespace casadi
