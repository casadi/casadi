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


#include "lace_modelica.hpp"
#include <lacemodelica/lacemodelica.h>

namespace casadi {

extern "C"
int CASADI_MODELICAPARSER_LACEMODELICA_EXPORT
casadi_register_modelicaparser_lacemodelica(ModelicaParserInternal::Plugin* plugin) {
  plugin->creator = LaceModelica::creator;
  plugin->name = "lacemodelica";
  plugin->doc = LaceModelica::meta_doc.c_str();
  plugin->version = CASADI_VERSION;
  return 0;
}

extern "C"
void CASADI_MODELICAPARSER_LACEMODELICA_EXPORT casadi_load_modelicaparser_lacemodelica() {
  ModelicaParserInternal::registerPlugin(casadi_register_modelicaparser_lacemodelica);
}

LaceModelica::LaceModelica()  {
}

LaceModelica::~LaceModelica() {
}

void LaceModelica::parse(const std::string& filename, const std::string& output_dir) {
  lacemodelica_status_t status = lacemodelica_process_bmo(filename.c_str(), output_dir.c_str());
  casadi_assert(status == LACEMODELICA_SUCCESS,
    "lacemodelica_process_bmo failed: " + std::string(lacemodelica_status_string(status)));
}

} // namespace casadi
