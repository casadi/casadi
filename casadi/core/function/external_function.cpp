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


#include "external_function_internal.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
//#include <dlfcn.h>

namespace casadi {

using namespace std;

  ExternalFunction::ExternalFunction() {
  }

  ExternalFunction::ExternalFunction(const std::string& bin_name) {
    assignNode(new ExternalFunctionInternal(bin_name));
  }

  ExternalFunctionInternal* ExternalFunction::operator->() {
    return static_cast<ExternalFunctionInternal*>(Function::operator->());
  }

  const ExternalFunctionInternal* ExternalFunction::operator->() const {
    return static_cast<const ExternalFunctionInternal*>(Function::operator->());
  }

  bool ExternalFunction::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const ExternalFunctionInternal*>(ptr)!=0;
  }

} // namespace casadi
