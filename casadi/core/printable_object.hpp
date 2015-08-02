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


#ifndef CASADI_PRINTABLE_OBJECT_HPP
#define CASADI_PRINTABLE_OBJECT_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <streambuf>
#include <vector>

#include "casadi_common.hpp"
#include "casadi_types.hpp"

namespace casadi {

  /** \brief Base class for objects that have a natural string representation
      \author Joel Andersson
      \date 2010-2014
  */
  template<class Derived>
  class CASADI_EXPORT PrintableObject {
  public:

    /// Return a string with a description (for SWIG)
    std::string getDescription() const {
      std::stringstream ss;
      static_cast<const Derived*>(this)->print(ss, false);
      return ss.str();
    }

    /// Return a string with a representation (for SWIG)
    std::string getRepresentation() const {
      std::stringstream ss;
      static_cast<const Derived*>(this)->repr(ss, false);
      return ss.str();
    }

#ifndef SWIG
    /// Print a representation of the object to a stream (shorthand)
    inline friend
      std::ostream& operator<<(std::ostream &stream, const PrintableObject<Derived>& obj) {
      static_cast<const Derived&>(obj).repr(stream, false);
      return stream;
    }
#endif // SWIG
/**
\ingroup expression_tools
@{
*/

#if !defined(SWIG) || defined(DOXYGEN)
    /// Return a string with a description of the object, cf. str(Object) in Python
    inline friend std::string str(const PrintableObject<Derived>& obj) {
      return obj.getDescription();
    }

    /// Return a string with a representation of the object, cf. repr(Object) in Python
    inline friend std::string repr(const PrintableObject<Derived>& obj) {
      return obj.getRepresentation();
    }
/** @} */
#endif //

  };
} // namespace casadi


#endif // CASADI_PRINTABLE_OBJECT_HPP
