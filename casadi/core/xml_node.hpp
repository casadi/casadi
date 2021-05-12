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


#ifndef CASADI_XML_NODE_HPP
#define CASADI_XML_NODE_HPP

#include <string>
#include <vector>
#include <map>
#include "exception.hpp"
#include "casadi_common.hpp"

/// \cond INTERNAL

namespace casadi {

struct CASADI_EXPORT XmlNode {
  // All attributes
  std::map<std::string, std::string> attributes;

  // All children
  std::vector<XmlNode> children;

  // Child index by name
  std::map<std::string, casadi_int> child_indices_;

  // Name of the node
  std::string name;

  // Comment
  std::string comment;

  // Text
  std::string text;

  /** \brief  Check if an attribute is present */
  bool has_attribute(const std::string& att_name) const;

  /** \brief  Add an attribute */
  void set_attribute(const std::string& att_name, const std::string& att);

  /** \brief  Get an attribute by its name */
  template<typename T>
  T attribute(const std::string& att_name) const {
    // Find the attribute, if any
    auto it = this->attributes.find(att_name);
    casadi_assert(it != this->attributes.end(), "Could not find attribute " + att_name);
    // Attribute found, read it
    T ret;
    read(it->second, &ret);
    return ret;
  }

  /** \brief  Get an attribute by its name, default value if not found */
  template<typename T>
  T attribute(const std::string& att_name, const T& def_att) const {
    // Find the attribute, if any
    auto it = this->attributes.find(att_name);
    if (it == this->attributes.end()) {
      // No such attribute, return default value
      return def_att;
    } else {
      // Attribute found, read it
      T ret;
      read(it->second, &ret);
      return ret;
    }
  }

  /** \brief  Get a reference to a child by its index */
  const XmlNode& operator[](size_t i) const { return this->children.at(i);}

  /** \brief  Get a reference to a child by its index */
  XmlNode& operator[](size_t i) { return this->children.at(i);}

  /** \brief  Get a reference to a child by its name */
  const XmlNode& operator[](const std::string& childname) const;

  /** \brief  Get a reference to a child by its name */
  XmlNode& operator[](const std::string& childname);

  /** \brief  Check if a child is present */
  bool has_child(const std::string& childname) const;

  /** \brief  Get the number of children */
  size_t size() const { return this->children.size();}

  /** \brief  Get value of text field */
  template<typename T>
    void get(T* val) const { read(this->text, val);}

  /** \brief  Read the string value of a string (i.e. copy) */
  static void read(const std::string& str, std::string* val);

  /** \brief  Read the boolean value of a string */
  static void read(const std::string& str, bool* val);

  /** \brief  Read the integer value of a string */
  static void read(const std::string& str, casadi_int* val);

  /** \brief  Read the double value of a string */
  static void read(const std::string& str, double* val);

  /** \brief  Read a vector of integer values of a string */
  static void read(const std::string& str, std::vector<casadi_int>* val);

  /** \brief Print to stream */
  CASADI_EXPORT friend std::ostream& operator<<(std::ostream &stream, const XmlNode& node);

  /** \brief  Dump representation */
  void dump(std::ostream &stream, casadi_int indent = 0) const;
};

} // namespace casadi
/// \endcond

#endif // CASADI_XML_NODE_HPP
