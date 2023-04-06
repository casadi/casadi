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


#include "xml_node.hpp"
#include "casadi_misc.hpp"

namespace casadi {

bool XmlNode::has_attribute(const std::string& att_name) const {
  return this->attributes.find(att_name) != this->attributes.end();
}

bool XmlNode::has_child(const std::string& childname) const {
  // Linear search over the elements
  for (auto& c : this->children) {
    if (c.name == childname) return true;
  }
  // Not found
  return false;
}

XmlNode& XmlNode::operator[](const std::string& childname) {
  // Linear search over the elements
  auto it = this->children.begin();
  for (; it != this->children.end(); ++it) {
    if (it->name == childname) break;
  }
  // Check that the child was indeed found
  casadi_assert(it != this->children.end(), "Could not find " + childname);
  // Return a reference to the child
  return *it;
}

const XmlNode& XmlNode::operator[](const std::string& childname) const {
  return const_cast<XmlNode*>(this)->operator[](childname); // NOLINT
}

void XmlNode::set_attribute(const std::string& attribute_name, const std::string& attribute) {
  this->attributes[attribute_name] = attribute;
}

void XmlNode::set_attribute(const std::string& att_name, const std::vector<casadi_int>& att) {
  std::stringstream ss;
  if (!att.empty()) {
    ss << att.at(0);
    for (size_t i = 1; i < att.size(); ++i) ss << " " << att.at(i);
  }
  return set_attribute(att_name, ss.str());
}

void XmlNode::set_attribute(const std::string& att_name, double att) {
  std::stringstream ss;
  ss << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1) << att;
  set_attribute(att_name, ss.str());
}

std::ostream& operator<<(std::ostream &stream, const XmlNode& node) {
  node.dump(stream);
  return stream;
}

void XmlNode::dump(std::ostream &stream, casadi_int indent) const {
  // Print name
  stream << std::string(indent, ' ') << "Node: " << this->name << std::endl;

  // Print comment
  if (!this->comment.empty()) {
    stream << std::string(indent, ' ') << "----- comment starts ----- "  << std::endl;
    stream << this->comment << std::endl;
    stream << std::string(indent, ' ') << "----- comment ends ----- "  << std::endl;
  }

  // Print text
  if (!this->text.empty())
    stream << std::string(indent+2, ' ') << "Text: " << this->text << std::endl;

  // Print attributes
  for (auto it = this->attributes.begin(); it != this->attributes.end(); ++it)
    stream << std::string(indent+2, ' ') << "attribute " << it->first
      << " = " << it->second << std::endl;

  // Print Children
  for (casadi_int i=0; i < size(); ++i) {
    stream << std::string(indent, ' ') << "Child " << i << ":" << std::endl;
    (*this)[i].dump(stream, indent+2);
  }
}

void XmlNode::read(const std::string& str, std::string* val) {
  *val = str;
}

void XmlNode::read(const std::string& str, bool* val) {
  if (str == "true") {
    *val = true;
  } else if (str == "false") {
    *val = false;
  } else {
    casadi_error("XML argument not 'true' or 'false'");
  }
}

void XmlNode::read(const std::string& str, casadi_int* val) {
  std::istringstream buffer(str);
  buffer >> *val;
}

void XmlNode::read(const std::string& str, double* val) {
  std::istringstream buffer(str);
  buffer >> *val;
}

void XmlNode::read(const std::string& str, std::vector<casadi_int>* val) {
  val->clear();
  std::istringstream buffer(str);
  while (true) {
    casadi_int v;
    buffer >> v;
    if (buffer.fail()) break;
    val->push_back(v);
  }
}

void XmlNode::read(const std::string& str, std::vector<std::string>* val) {
  val->clear();
  std::istringstream buffer(str);
  while (true) {
    std::string v;
    buffer >> v;
    if (buffer.fail()) break;
    val->push_back(v);
  }
}

std::vector<std::string> XmlNode::child_names() const {
  std::vector<std::string> ret;
  ret.reserve(this->children.size());
  for (auto& c : this->children) ret.push_back(c.name);
  return ret;
}

std::vector<std::string> XmlNode::attribute_names() const {
  std::vector<std::string> ret;
  ret.reserve(this->attributes.size());
  for (auto& a : this->attributes) ret.push_back(a.first);
  return ret;
}

} // namespace casadi
