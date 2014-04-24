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

#include "xml_node.hpp"
#include "external_packages/tinyxml/tinyxml.h"
#include "casadi/core/casadi_exception.hpp"

using namespace std;
namespace casadi {

XMLNode::XMLNode() {
}

XMLNode::~XMLNode() {
}

bool XMLNode::hasAttribute(const string& attribute_name) const {
  map<string, string>::const_iterator it = attributes_.find(attribute_name);
  return it!=attributes_.end();
}

XMLNode& XMLNode::operator[](int i) {
  casadi_assert_message(i>=0 && i < size(),
                        "XMLNode::operator[]: index out of bounds for element " << i
                        << " of node " << getName());

  return children_.at(i);
}

const XMLNode& XMLNode::operator[](int i) const {
  return const_cast<XMLNode*>(this)->operator[](i);
}

bool XMLNode::hasChild(const string& childname) const {
  map<string, int>::const_iterator it = child_indices_.find(childname);
  return it!=child_indices_.end();
}

XMLNode& XMLNode::operator[](const string& childname) {
  // Find the child
  map<string, int>::const_iterator it = child_indices_.find(childname);

  // check that the child was indeed found
  if (it == child_indices_.end()) {
    throw CasadiException("Error in XMLNode::operator[]: could not find " + childname);
  }

  // Return an index to the child
  return children_[it->second];
}

const XMLNode& XMLNode::operator[](const string& childname) const {
  return const_cast<XMLNode*>(this)->operator[](childname);
}

void XMLNode::setAttribute(const string& attribute_name, const string& attribute) {
  attributes_[attribute_name] = attribute;
}

void XMLNode::addNode(TiXmlNode* n) {
  if (!n) throw CasadiException("Error in XMLNode::addNode: Node is 0");

  // Save name
  setName(n->Value());

  // Save attributes
  int type = n->Type();
  if (type == TiXmlNode::ELEMENT) {
    if (n->ToElement()!=0) {
      for (TiXmlAttribute* pAttrib=n->ToElement()->FirstAttribute();
          pAttrib;
          pAttrib=pAttrib->Next()) {
          setAttribute(pAttrib->Name(), pAttrib->Value());
      }
    }
  } else if (type == TiXmlNode::DOCUMENT) {
    // do nothing
  } else {
    throw CasadiException("XMLNode::addNode");
  }

  // Count the number of children
  int num_children = 0;
  for (TiXmlNode* child = n->FirstChild(); child != 0; child= child->NextSibling()) {
    num_children++;
  }
  children_.reserve(num_children);

  // add children
  int ch = 0;
  for (TiXmlNode* child = n->FirstChild(); child != 0; child= child->NextSibling(), ++ch) {
      int childtype = child->Type();

      if (childtype == TiXmlNode::ELEMENT) {
        XMLNode newnode;
        newnode.addNode(child);
        children_.push_back(newnode);
        child_indices_[newnode.getName()] = ch;
      } else if (childtype == TiXmlNode::COMMENT) {
        comment_ = child->Value();
      } else if (childtype == TiXmlNode::TEXT) {
        text_ = child->ToText()->Value();
      } else if (childtype == TiXmlNode::DECLARATION) {
        cout << "Warning: Skipped TiXmlNode::DECLARATION" << endl;
      } else {
        throw CasadiException("Error in XMLNode::addNode: Unknown node type");
      }
  }
}

ostream& operator<<(ostream &stream, const XMLNode& node) {
  node.dump(stream);
  return stream;
}

int XMLNode::size() const {
  return children_.size();
}

const string& XMLNode::getName() const {
  return name_;
}

void XMLNode::setName(const string& name) {
  name_ = name;
}

void XMLNode::dump(ostream &stream, int indent) const {
  // Print name
  stream << string(indent, ' ') << "Node: " << name_ << endl;

  // Print comment
  if (!comment_.empty()) {
    stream << string(indent, ' ') << "----- comment starts ----- "  << endl;
    stream << comment_ << endl;
    stream << string(indent, ' ') << "----- comment ends ----- "  << endl;
  }

  // Print text
  if (!text_.empty())
    stream << string(indent+2, ' ') << "Text: " << text_ << endl;

  // Print attributes
  for (map<string, string>::const_iterator it=attributes_.begin(); it != attributes_.end(); ++it)
    stream << string(indent+2, ' ') << "attribute " << it->first << " = " << it->second << endl;

  // Print Children
  for (int i=0; i<size(); ++i) {
    stream << string(indent, ' ') << "Child " << i << ":" << endl;
    (*this)[i].dump(stream, indent+2);
  }
}

bool XMLNode::checkName(const string& str) const {
  return name_.compare(str) == 0;
}

void XMLNode::readString(const std::string& str, std::string& val) {
  val = str;
}

void XMLNode::readString(const std::string& str, bool& val) {
  if (str.compare("true")==0)
    val = true;
  else if (str.compare("false")==0)
    val = false;
  else
    throw CasadiException("XML argument not true or false");
}

void XMLNode::readString(const std::string& str, int& val) {
  std::istringstream buffer(str);
  buffer >> val;
}

void XMLNode::readString(const std::string& str, double& val) {
  std::istringstream buffer(str);
  buffer >> val;
}

} // namespace casadi
