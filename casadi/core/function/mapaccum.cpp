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


#include "mapaccum_internal.hpp"

namespace casadi {
  using namespace std;

  MapAccum::MapAccum() {
  }

  MapAccum::MapAccum(const std::string& name, const Function& f,
                 int n,
                 const std::vector<bool>& input_accum,
                 const std::vector<int>& output_accum,
                 bool reverse,
                 const Dict& opts) {
    assignNode(new MapAccumInternal(name, f, n, input_accum, output_accum, reverse));
    setOption(opts);
    init();
  }

  MapAccumInternal* MapAccum::operator->() {
    return static_cast<MapAccumInternal*>(Function::operator->());
  }

  const MapAccumInternal* MapAccum::operator->() const {
    return static_cast<const MapAccumInternal*>(Function::operator->());
  }

  bool MapAccum::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const MapAccumInternal*>(ptr)!=0;
  }

} // namespace casadi
