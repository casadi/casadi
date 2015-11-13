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


#include "simulator_internal.hpp"
#include "ivpsol.hpp"
#include "../std_vector_tools.hpp"

using namespace std;
namespace casadi {


  SimulatorInternal::SimulatorInternal(const std::string& name, const Function& integrator)
    : FunctionInternal(name), integrator_(integrator) {
    grid_ = dynamic_cast<Ivpsol*>(integrator_.get())->grid_;
  }

  SimulatorInternal::~SimulatorInternal() {
  }

  void SimulatorInternal::init() {
    ischeme_ = Function::ivpsol_in();
    oscheme_ = Function::ivpsol_out();

    // Call base class method
    FunctionInternal::init();
    alloc(integrator_);
  }

  void SimulatorInternal::eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    integrator_(arg, res, iw, w, mem);
  }

} // namespace casadi
