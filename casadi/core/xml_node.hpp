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

  // Name of the node
  std::string name;

  // Comment
  std::string comment;

  // Line number
  casadi_int line;

  // Text
  std::string text;

  /** \brief  Check if an attribute is present

      \identifier{vc} */
  bool has_attribute(const std::string& att_name) const;

  /** \brief  Add an attribute

      \identifier{vd} */
  void set_attribute(const std::string& att_name, const std::string& att);

  /** \brief Add an integer attribute

      \identifier{252} */
  void set_attribute(const std::string& att_name, casadi_int att) {
    set_attribute(att_name, std::to_string(att));
  }

  /** \brief Add an integer attribute

      \identifier{253} */
  void set_attribute(const std::string& att_name, double att);

  /** \brief Add a vector attribute

      \identifier{254} */
  void set_attribute(const std::string& att_name, const std::vector<casadi_int>& att);

  /** \brief  Names of children

      \identifier{ve} */
  std::vector<std::string> child_names() const;

  /** \brief  Names of attributes

      \identifier{vf} */
  std::vector<std::string> attribute_names() const;

  /** \brief  Get an attribute by its name

      \identifier{vg} */
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

  /** \brief  Get an attribute by its name, default value if not found

      \identifier{vh} */
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

  /** \brief  Get a reference to a child by its index

      \identifier{vi} */
  const XmlNode& operator[](size_t i) const { return this->children.at(i);}

  /** \brief  Get a reference to a child by its index

      \identifier{vj} */
  XmlNode& operator[](size_t i) { return this->children.at(i);}

  /** \brief  Get a reference to a child by its name

      \identifier{vk} */
  const XmlNode& operator[](const std::string& childname) const;

  /** \brief  Get a reference to a child by its name

      \identifier{vl} */
  XmlNode& operator[](const std::string& childname);

  /** \brief  Check if a child is present

      \identifier{vm} */
  bool has_child(const std::string& childname) const;

  /** \brief  Get the number of children

      \identifier{vn} */
  size_t size() const { return this->children.size();}

  /** \brief  Get value of text field

      \identifier{vo} */
  template<typename T>
  void get(T* val) const { read(this->text, val);}

  /** \brief  Read the string value of a string (i.e. copy)

      \identifier{vp} */
  static void read(const std::string& str, std::string* val);

  /** \brief  Read the boolean value of a string

      \identifier{vq} */
  static void read(const std::string& str, bool* val);

  /** \brief  Read the integer value of a string

      \identifier{vr} */
  static void read(const std::string& str, casadi_int* val);

  /** \brief  Read the size_t value of a string

      \identifier{29u} */
  static void read(const std::string& str, size_t* val);

  /** \brief  Read the double value of a string

      \identifier{vs} */
  static void read(const std::string& str, double* val);

  /** \brief  Read a vector of integer values of a string

      \identifier{vt} */
  static void read(const std::string& str, std::vector<casadi_int>* val);

  /** \brief  Read a vector of string values of a string

      \identifier{vu} */
  static void read(const std::string& str, std::vector<std::string>* val);

  /** \brief Print to stream

      \identifier{vv} */
  CASADI_EXPORT friend std::ostream& operator<<(std::ostream &stream, const XmlNode& node);

  /** \brief  Dump representation

      \identifier{vw} */
  void dump(std::ostream &stream, casadi_int indent = 0) const;
};

} // namespace casadi
/// \endcond

#endif // CASADI_XML_NODE_HPP
