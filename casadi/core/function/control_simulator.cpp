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


#include "control_simulator.hpp"
#include "control_simulator_internal.hpp"

using namespace std;
namespace casadi {

  ControlSimulator::ControlSimulator() {
  }

  ControlSimulator::ControlSimulator(const Function& dae, const Function& output_fcn,
                                     const vector<double>& grid) {
    assignNode(new ControlSimulatorInternal(dae, output_fcn, grid));
  }

  ControlSimulator::ControlSimulator(const Function& dae, const Function& output_fcn,
                                     const Matrix<double>& grid) {
    casadi_assert_message(grid.isVector(), "ControlSimulator::ControlSimulator: grid must be a "
                          "column vector, but got " << grid.dimString());
    casadi_assert_message(grid.isDense(), "ControlSimulator::ControlSimulator: grid must be dense, "
                          "but got " << grid.dimString());
    assignNode(new ControlSimulatorInternal(dae, output_fcn, grid.data()));
  }

  ControlSimulator::ControlSimulator(const Function& dae, const vector<double>& grid) {
    assignNode(new ControlSimulatorInternal(dae, Function(), grid));
  }

  ControlSimulator::ControlSimulator(const Function& dae, const Matrix<double>& grid) {
    casadi_assert_message(grid.isVector(), "ControlSimulator::ControlSimulator: grid must be a "
                          "column vector, but got " << grid.dimString());
    casadi_assert_message(grid.isDense(), "ControlSimulator::ControlSimulator: grid must be dense, "
                          "but got " << grid.dimString());
    assignNode(new ControlSimulatorInternal(dae, Function(), grid.data()));
  }

  ControlSimulatorInternal* ControlSimulator::operator->() {
    return static_cast<ControlSimulatorInternal*>(Function::operator->());
  }

  const ControlSimulatorInternal* ControlSimulator::operator->() const {
    return static_cast<const ControlSimulatorInternal*>(Function::operator->());
  }

  bool ControlSimulator::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const ControlSimulatorInternal*>(ptr)!=0;
  }

  std::vector<double> ControlSimulator::getMinorT() const {
    casadi_assert(!isNull());
    return dynamic_cast<const ControlSimulatorInternal*>(get())->grid_;
  }

  Matrix<double> ControlSimulator::getMinorU() const {
    return dynamic_cast<const ControlSimulatorInternal*>(get())->getVFine();
  }

  std::vector<int> ControlSimulator::getMajorIndex() const {
    return dynamic_cast<const ControlSimulatorInternal*>(get())->getCoarseIndex();
  }

} // namespace casadi

