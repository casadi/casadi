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


#ifndef CASADI_ENUM_HPP
#define CASADI_ENUM_HPP

#include <string>
#include <sstream>

#include "exception.hpp"

#include <casadi/core/casadi_export.h>

/// \cond INTERNAL
namespace casadi {

/// Helper class: Specify number of entries in an enum
template<typename T>
struct enum_traits {
  static const size_t n_enum = static_cast<size_t>(T::NUMEL);
};

/// Helper function: Check if enum exists
template<typename T>
bool has_enum(const std::string& s) {
  // Look for a match
  for (size_t i = 0; i < enum_traits<T>::n_enum; ++i) {
    if (s == to_string(static_cast<T>(i))) return true;
  }
  // Not found
  return false;
}

/// Helper function: Convert string to enum
template<typename T>
T to_enum(const std::string& s, const std::string& s_def = "") {
  // Default value, if empty string
  if (s.empty() && !s_def.empty()) return to_enum<T>(s_def);
  // Linear search over permitted values
  for (size_t i = 0; i < enum_traits<T>::n_enum; ++i) {
    if (s == to_string(static_cast<T>(i))) return static_cast<T>(i);
  }
  // Informative error message
  std::stringstream ss;
  ss << "No such enum: '" << s << "'. Permitted values: ";
  for (size_t i = 0; i < enum_traits<T>::n_enum; ++i) {
    // Separate strings
    if (i > 0) ss << ", ";
    // Print enum name
    ss << "'" << to_string(static_cast<T>(i)) << "'";
  }
  casadi_error(ss.str());
  return static_cast<T>(enum_traits<T>::n_enum);  // never reached
}

/// Helper function: Get all fields
template<typename T>
std::vector<std::string> enum_names() {
  std::vector<std::string> r(enum_traits<T>::n_enum);
  for (size_t i = 0; i < enum_traits<T>::n_enum; ++i)
    r[i] = to_string(static_cast<T>(i));
  return r;
}

} // namespace casadi
/// \endcond

#endif // CASADI_ENUM_HPP
