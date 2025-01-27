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


#ifndef CASADI_PRINTABLE_HPP
#define CASADI_PRINTABLE_HPP

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnonnull"
#endif
#include <iostream>
#if __GNUC__
#pragma GCC diagnostic pop
#endif

#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

#include "casadi_common.hpp"

namespace casadi {

  /** \brief Empty Base

      This class is extended in SWIG.

      \identifier{7h} */
  struct CASADI_EXPORT PrintableCommon {};

#ifndef SWIG
  /** \brief Base class for objects that have a natural string representation

      \author Joel Andersson
      \date 2010-2014

      \identifier{7i} */
  template<class Derived>
  class CASADI_EXPORT Printable : public PrintableCommon {
  public:
    /// Print a string representation of the object to a stream
    inline friend
      std::ostream& operator<<(std::ostream &stream, const Derived& obj) {
      obj.disp(stream);
      return stream;
    }

    /// Get string representation with type information
    inline friend std::string repr(const Derived& obj) {
      return Derived::type_name() + "(" + obj.get_str() + ")";
    }
  };
#endif // SWIG

} // namespace casadi


#endif // CASADI_PRINTABLE_HPP
