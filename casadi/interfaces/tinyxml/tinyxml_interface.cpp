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

#include "tinyxml_interface.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_XMLFILE_TINYXML_EXPORT
  casadi_register_xmlfile_tinyxml(XmlFileInternal::Plugin* plugin) {
    plugin->creator = TinyXmlInterface::creator;
    plugin->name = "tinyxml";
    plugin->doc = TinyXmlInterface::meta_doc.c_str();
    plugin->version = 20;
    return 0;
  }

  extern "C"
  void CASADI_XMLFILE_TINYXML_EXPORT casadi_load_xmlfile_tinyxml() {
    XmlFileInternal::registerPlugin(casadi_register_xmlfile_tinyxml);
  }

  TinyXmlInterface::TinyXmlInterface() : XmlFileInternal() {
  }

  TinyXmlInterface::~TinyXmlInterface() {
  }

  TinyXmlInterface* TinyXmlInterface::clone() const {
    return new TinyXmlInterface();
  }

} // namespace casadi
