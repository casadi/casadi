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


#include "xml_node.hpp"
#include "casadi_misc.hpp"

namespace casadi {

bool XmlNode::has_attribute(const std::string& att_name) const {
  return attributes_.find(att_name) != attributes_.end();
}

bool XmlNode::has_child(const std::string& childname) const {
  auto it = child_indices_.find(childname);
  return it!=child_indices_.end();
}

XmlNode& XmlNode::operator[](const std::string& childname) {
  // Find the child
  auto it = child_indices_.find(childname);

  // check that the child was indeed found
  if (it == child_indices_.end()) {
    casadi_error("could not find " + childname);
  }

  // Return an index to the child
  return children_[it->second];
}

const XmlNode& XmlNode::operator[](const std::string& childname) const {
  return const_cast<XmlNode*>(this)->operator[](childname); // NOLINT
}

void XmlNode::set_attribute(const std::string& attribute_name, const std::string& attribute) {
  attributes_[attribute_name] = attribute;
}

std::ostream& operator<<(std::ostream &stream, const XmlNode& node) {
  node.dump(stream);
  return stream;
}


void XmlNode::dump(std::ostream &stream, casadi_int indent) const {
  // Print name
  stream << std::string(indent, ' ') << "Node: " << this->name << std::endl;

  // Print comment
  if (!comment_.empty()) {
    stream << std::string(indent, ' ') << "----- comment starts ----- "  << std::endl;
    stream << comment_ << std::endl;
    stream << std::string(indent, ' ') << "----- comment ends ----- "  << std::endl;
  }

  // Print text
  if (!text_.empty())
    stream << std::string(indent+2, ' ') << "Text: " << text_ << std::endl;

  // Print attributes
  for (auto it=attributes_.begin(); it != attributes_.end(); ++it)
    stream << std::string(indent+2, ' ') << "attribute " << it->first
      << " = " << it->second << std::endl;

  // Print Children
  for (casadi_int i=0; i < size(); ++i) {
    stream << std::string(indent, ' ') << "Child " << i << ":" << std::endl;
    (*this)[i].dump(stream, indent+2);
  }
}

void XmlNode::readString(const std::string& str, std::string& val) {
  val = str;
}

void XmlNode::readString(const std::string& str, bool& val) {
  if (str=="true")
    val = true;
  else if (str=="false")
    val = false;
  else
    throw CasadiException("XML argument not true or false");
}

void XmlNode::readString(const std::string& str, casadi_int& val) {
  std::istringstream buffer(str);
  buffer >> val;
}

void XmlNode::readString(const std::string& str, double& val) {
  std::istringstream buffer(str);
  buffer >> val;
}

void XmlNode::readString(const std::string& str, std::vector<casadi_int>& val) {
  val.clear();
  std::istringstream buffer(str);
  while (true) {
    casadi_int v;
    buffer >> v;
    if (buffer.fail()) break;
    val.push_back(v);
  }
}

} // namespace casadi
