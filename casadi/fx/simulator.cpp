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

#include "simulator.hpp"
#include "simulator_internal.hpp"
#include "integrator_internal.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"

using namespace std;
namespace CasADi{

Simulator::Simulator(){
}

Simulator::Simulator(const Integrator& integrator, const FX& output_fcn, const vector<double>& grid){
  assignNode(new SimulatorInternal(integrator,output_fcn,grid));
}

Simulator::Simulator(const Integrator& integrator, const vector<double>& grid){
  // Create a dummy function (returns the whole state)
  SXMatrix t = symbolic("t");
  SXMatrix x = symbolic("x",integrator->nx_);
  SXMatrix p = symbolic("p",integrator->np_);
  vector<SXMatrix> arg(OUTPUT_NUM_IN);  arg[OUTPUT_T] = t; arg[OUTPUT_X] = x; arg[OUTPUT_P] = p;
  
  // Create the output function
  SXFunction output_fcn(arg,vector<SXMatrix>(1,x));
  assignNode(new SimulatorInternal(integrator,output_fcn,grid));
}

SimulatorInternal* Simulator::operator->(){
  return (SimulatorInternal*)(FX::operator->());
}

const SimulatorInternal* Simulator::operator->() const{
   return (const SimulatorInternal*)(FX::operator->()); 
}

bool Simulator::checkNode() const{
  return dynamic_cast<const SimulatorInternal*>(get())!=0;
}


} // namespace CasADi

