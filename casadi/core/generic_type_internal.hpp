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


#ifndef CASADI_GENERIC_TYPE_INTERNAL_HPP
#define CASADI_GENERIC_TYPE_INTERNAL_HPP

#include "generic_type.hpp"
#include "std_vector_tools.hpp"

/// \cond INTERNAL

namespace casadi {

  template<typename T>
  class CASADI_CORE_EXPORT GenericTypeInternal : public SharedObjectNode {
    public:
      explicit GenericTypeInternal(const T& d) : d_(d) {}
      virtual ~GenericTypeInternal() {}
      virtual GenericTypeInternal<T>* clone() const { return new GenericTypeInternal(d_);}
      virtual void print(std::ostream &stream) const { stream << d_; }
      T d_;
  };

  //template<> GenericTypeInternal<bool>::type = OT_BOOLEAN;

  // Implementations of public functions
  template<typename T>
  bool GenericType::is_a() const {
    return dynamic_cast<const GenericTypeInternal<T>*>(get()) != 0;
  }

} // namespace casadi
/// \endcond

#endif // CASADI_GENERIC_TYPE_INTERNAL_HPP
