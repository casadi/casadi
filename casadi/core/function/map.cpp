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


#include "map_internal.hpp"


namespace casadi {
  using namespace std;

  Map::Map() {
  }

  Map::Map(const std::string& name, const Function& f, int n, const Dict& opts) {
    assignNode(MapBase::create(f, n, opts));
    setOption("name", name);
    setOption(opts);
    init();
  }

  Map::Map(const std::string& name, const Function& f,
                 int n,
                 const std::vector<bool> &repeat_in,
                 const std::vector<bool> &repeat_out,
                 const Dict& opts) {
    Dict options = opts;
    std::vector<int> reduced_inputs;
    std::vector<int> reduced_outputs;

    for (int i=0;i<repeat_in.size();++i) {
      if (!repeat_in[i]) reduced_inputs.push_back(i);
    }
    for (int i=0;i<repeat_out.size();++i) {
      if (!repeat_out[i]) reduced_outputs.push_back(i);
    }

    options["reduced_inputs"] = reduced_inputs;
    options["reduced_outputs"] = reduced_outputs;

    assignNode(MapBase::create(f, n, options));
    setOption("name", name);
    setOption(opts);
    init();
  }

  MapBase* Map::operator->() {
    return static_cast<MapBase*>(Function::operator->());
  }

  const MapBase* Map::operator->() const {
    return static_cast<const MapBase*>(Function::operator->());
  }

  bool Map::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const MapBase*>(ptr)!=0;
  }

} // namespace casadi
