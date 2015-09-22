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


#include "switch_internal.hpp"

namespace casadi {
  using namespace std;

  Switch::Switch() {
  }

  Switch::Switch(const std::string& name, const std::vector<Function>& f,
                 const Function& f_def, const Dict& opts) {
    assignNode(new SwitchInternal(f, f_def));
    setOption("name", name);
    setOption(opts);
    init();
  }

  Switch::Switch(const std::string& name, const Function& f_true,
                 const Function& f_false, const Dict& opts) {
    assignNode(new SwitchInternal(vector<Function>(1, f_false), f_true));
    setOption("name", name);
    setOption(opts);
    init();
  }

  SwitchInternal* Switch::operator->() {
    return static_cast<SwitchInternal*>(Function::operator->());
  }

  const SwitchInternal* Switch::operator->() const {
    return static_cast<const SwitchInternal*>(Function::operator->());
  }

  bool Switch::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SwitchInternal*>(ptr)!=0;
  }

} // namespace casadi
