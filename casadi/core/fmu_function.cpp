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


#include "fmu_function_impl.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

namespace casadi {

Function fmu_function(const std::string& name, const std::string& guid,
    const std::string& resource_loc,
    const std::vector<std::vector<casadi_int>>& id_in,
    const std::vector<std::vector<casadi_int>>& id_out,
    const Dict& opts) {
  return Function::create(new FmuFunction(name, guid, resource_loc, id_in, id_out), opts);
}

FmuFunction::FmuFunction(const std::string& name, const std::string& guid,
    const std::string& resource_loc,
    const std::vector<std::vector<casadi_int>>& id_in,
    const std::vector<std::vector<casadi_int>>& id_out)
  : FunctionInternal(name),
    guid_(guid), resource_loc_(resource_loc), id_in_(id_in), id_out_(id_out) {
}

FmuFunction::~FmuFunction() {
  clear_mem();
}

void FmuFunction::init(const Dict& opts) {
  // Call the initialization method of the base class
  FunctionInternal::init(opts);
}

} // namespace casadi
