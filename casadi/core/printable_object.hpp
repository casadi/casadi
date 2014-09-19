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

#ifndef CASADI_PRINTABLE_OBJECT_HPP
#define CASADI_PRINTABLE_OBJECT_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <streambuf>
#include <vector>

#include "casadi_common.hpp"

namespace casadi {

  /** \brief Base class for objects that have a natural string representation
      \author Joel Andersson
      \date 2010
  */
  class CASADI_CORE_EXPORT PrintableObject {
  public:

#ifndef SWIG
    /// Print a description of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Print a representation of the object to a stream (shorthand)
    CASADI_CORE_EXPORT friend std::ostream& operator<<(std::ostream &stream,
                                                           const PrintableObject& obj);

    /// Return a string with a description of the object, cf. str(Object) in Python
    CASADI_CORE_EXPORT friend std::string str(const PrintableObject& obj);

    /// Return a string with a representation of the object, cf. repr(Object) in Python
    CASADI_CORE_EXPORT friend std::string repr(const PrintableObject& obj);

#endif // SWIG

    /// Return a string with a representation (for SWIG)
    std::string getRepresentation() const;

    /// Return a string with a description (for SWIG)
    std::string getDescription() const;
  };

#ifdef SWIG
  %extend PrintableObject {
    std::string __str__()  { return $self->getDescription(); }
    std::string __repr__()  { return $self->getRepresentation(); }
  }
#endif // SWIG

} // namespace casadi


#endif // CASADI_PRINTABLE_OBJECT_HPP

