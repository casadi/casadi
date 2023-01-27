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


#include "tinyxml_interface.hpp"

namespace casadi {

  extern "C"
  int CASADI_XMLFILE_TINYXML_EXPORT
  casadi_register_xmlfile_tinyxml(XmlFileInternal::Plugin* plugin) {
    plugin->creator = TinyXmlInterface::creator;
    plugin->name = "tinyxml";
    plugin->doc = TinyXmlInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    return 0;
  }

  extern "C"
  void CASADI_XMLFILE_TINYXML_EXPORT casadi_load_xmlfile_tinyxml() {
    XmlFileInternal::registerPlugin(casadi_register_xmlfile_tinyxml);
  }

  TinyXmlInterface::TinyXmlInterface()  {
  }

  TinyXmlInterface::~TinyXmlInterface() {
  }

  XmlNode TinyXmlInterface::parse(const std::string& filename) {
    XMLError err = doc_.LoadFile(filename.c_str());
    casadi_assert(!err, "Cound not open " + filename);
    return addNode(&doc_);
  }

  XmlNode TinyXmlInterface::addNode(TiXmlNode* n) {
    if (!n) casadi_error("Node is 0");
    XmlNode ret;

    ret.line = n->GetLineNum();

    // Save name
    if (n->Value()) {
      ret.name = n->Value();
    }

    // Save attributes
    if (n->ToElement()) {    
      for (const TiXmlAttribute* pAttrib=n->ToElement()->FirstAttribute();
                                 pAttrib; pAttrib=pAttrib->Next()) {
        ret.set_attribute(pAttrib->Name(), pAttrib->Value());
      }
    } else if (n->ToDocument()) {
      // do nothing
    } else {
      casadi_error("TinyXmlInterface::addNode");
    }

    // Count the number of children
    casadi_int num_children = 0;
    for (TiXmlNode* child = n->FirstChild(); child != nullptr; child = child->NextSibling()) {
      num_children++;
    }
    ret.children.reserve(num_children);

    // add children
    for (TiXmlNode* child = n->FirstChild(); child != nullptr; child = child->NextSibling()) {
      if (child->ToElement()) {
        ret.children.push_back(addNode(child));
      } else if (child->ToComment()) {
        ret.comment = child->Value();
      } else if (child->ToText()) {
        ret.text = child->ToText()->Value();
      } else if (child->ToDeclaration()) {
        // pass
      } else if (child->ToDocument()) {
        // pass
      } else {
        casadi_error("Unknown node type");
      }
    }

    // Note: Return value optimization
    return ret;
  }


} // namespace casadi
