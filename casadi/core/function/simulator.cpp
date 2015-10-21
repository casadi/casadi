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


#include "simulator.hpp"
#include "simulator_internal.hpp"

using namespace std;
namespace casadi {

  Simulator::Simulator() {
  }

  Simulator::Simulator(const std::string& name, const Function& integrator,
                       const Function& output_fcn,
                       const Matrix<double>& grid, const Dict& opts) {
    assignNode(new SimulatorInternal(name, integrator, output_fcn, grid));
    setOption(opts);
    init();
  }

  Simulator::Simulator(const std::string& name, const Function& integrator,
                       const Matrix<double>& grid, const Dict& opts) {
    assignNode(new SimulatorInternal(name, integrator, Function(), grid));
    setOption(opts);
    init();
  }

  SimulatorInternal* Simulator::operator->() {
    return static_cast<SimulatorInternal*>(Function::operator->());
  }

  const SimulatorInternal* Simulator::operator->() const {
    return static_cast<const SimulatorInternal*>(Function::operator->());
  }

  bool Simulator::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SimulatorInternal*>(ptr)!=0;
  }


} // namespace casadi
