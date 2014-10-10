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

  Simulator::Simulator(const Integrator& integrator, const Function& output_fcn,
                       const vector<double>& grid) {
    assignNode(new SimulatorInternal(integrator, output_fcn, grid));
  }

  Simulator::Simulator(const Integrator& integrator, const Function& output_fcn,
                       const Matrix<double>& grid) {
    casadi_assert_message(grid.isVector(),
                          "Simulator::Simulator: grid must be a column vector, but got "
                          << grid.dimString());
    casadi_assert_message(grid.isDense(),
                          "Simulator::Simulator: grid must be dense, but got "
                          << grid.dimString());
    assignNode(new SimulatorInternal(integrator, output_fcn, grid.data()));
  }

  Simulator::Simulator(const Integrator& integrator, const vector<double>& grid) {
    assignNode(new SimulatorInternal(integrator, Function(), grid));
  }

  Simulator::Simulator(const Integrator& integrator, const Matrix<double>& grid) {
    casadi_assert_message(grid.isVector(),
                          "Simulator::Simulator: grid must be column vector, but got "
                          << grid.dimString());
    casadi_assert_message(grid.isDense(),
                          "Simulator::Simulator: grid must be dense, but got "
                          << grid.dimString());
    assignNode(new SimulatorInternal(integrator, Function(), grid.data()));
  }


  SimulatorInternal* Simulator::operator->() {
    return static_cast<SimulatorInternal*>(Function::operator->());
  }

  const SimulatorInternal* Simulator::operator->() const {
    return static_cast<const SimulatorInternal*>(Function::operator->());
  }

  bool Simulator::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SimulatorInternal*>(ptr)!=0;
  }


} // namespace casadi
