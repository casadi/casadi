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

#include "shared_object.hpp"
#include "generic_shared_impl.hpp"

namespace casadi {


  std::string SharedObject::class_name() const {
    return (*this)->class_name();
  }

  void SharedObject::disp(std::ostream& stream, bool more) const {
    if (is_null()) {
      stream << "NULL";
    } else {
      (*this)->disp(stream, more);
    }
  }

  void SharedObject::print_ptr(std::ostream &stream) const {
    stream << get();
  }

  void WeakRefInternal::disp(std::ostream& stream, bool more) const {
    if (raw_==nullptr) {
      stream << "NULL";
    } else {
      raw_->disp(stream, more);
    }
  }

template class CASADI_EXPORT GenericShared< SharedObject, SharedObjectInternal >;
template class CASADI_EXPORT GenericWeakRef< SharedObject, SharedObjectInternal >;

} // namespace casadi
