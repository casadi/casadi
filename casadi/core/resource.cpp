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


#include "resource.hpp"
#include "resource_internal.hpp"

namespace casadi {

  Resource::Resource() {
  }

  Resource::Resource(const std::string& path) {
    // Pass through empty strings
    if (path.empty()) {
      own(new DirResource(path));
      return;
    }
    std::ifstream val(path);
    if (val) { // path is an existing file (hence not a directory)
      own(new ZipResource(path));
      return;
    } else { // path does not exist or is a directory
      own(new DirResource(path));
    }
  }

  Resource Resource::create(ResourceInternal *node) {
    Resource ret;
    ret.own(node);
    return ret;
  }

  const ResourceInternal* Resource::operator->() const {
    return static_cast<const ResourceInternal*>(SharedObject::operator->());
  }

  const ResourceInternal& Resource::operator*() const {
    return *static_cast<const ResourceInternal*>(get());
  }

  bool Resource::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const ResourceInternal*>(ptr)!=nullptr;
  }

  ResourceInternal* Resource::get() const {
    return static_cast<ResourceInternal*>(SharedObject::get());
  }

  const std::string& Resource::path() const {
    return (*this)->path();
  }

  void Resource::serialize(SerializingStream &s) const {
    return (*this)->serialize(s);
  }

  Resource Resource::deserialize(DeserializingStream& s) {
    return Resource::create(ResourceInternal::deserialize(s));
  }

} // namespace casadi
