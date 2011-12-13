/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "piecewise_simulator.hpp"
#include "piecewise_simulator_internal.hpp"

using namespace std;
namespace CasADi{

PiecewiseSimulator::PiecewiseSimulator(){
}

PiecewiseSimulator::PiecewiseSimulator(const FX& dae, const FX& output_fcn, const vector<double>& grid){
  assignNode(new PiecewiseSimulatorInternal(dae,output_fcn,grid));
}

PiecewiseSimulator::PiecewiseSimulator(const FX& dae, const vector<double>& grid){
  assignNode(new PiecewiseSimulatorInternal(dae,FX(),grid));
}

PiecewiseSimulatorInternal* PiecewiseSimulator::operator->(){
  return (PiecewiseSimulatorInternal*)(FX::operator->());
}

const PiecewiseSimulatorInternal* PiecewiseSimulator::operator->() const{
   return (const PiecewiseSimulatorInternal*)(FX::operator->()); 
}

bool PiecewiseSimulator::checkNode() const{
  return dynamic_cast<const PiecewiseSimulatorInternal*>(get())!=0;
}

std::vector<double> PiecewiseSimulator::getGrid() const { 
 	  casadi_assert(checkNode()); 
 	  return dynamic_cast<const PiecewiseSimulatorInternal*>(get())->grid_; 
} 

Matrix<double> PiecewiseSimulator::getVFine() const {
 	  return dynamic_cast<const PiecewiseSimulatorInternal*>(get())->getVFine(); 
}

} // namespace CasADi

