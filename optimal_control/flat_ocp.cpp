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

#include "flat_ocp_internal.hpp"

using namespace std;
namespace CasADi{
namespace OptimalControl{

FlatOCP::FlatOCP(){
}
    
FlatOCP::FlatOCP(const std::string& filename){
  assignNode(new FlatOCPInternal(filename));
}

void FlatOCP::parse(){
  (*this)->parse();
}

FlatOCPInternal* FlatOCP::operator->(){
  return (FlatOCPInternal*)(OptionsFunctionality::operator->());
}

const FlatOCPInternal* FlatOCP::operator->() const{
  return (const FlatOCPInternal*)(OptionsFunctionality::operator->());
}

bool FlatOCP::checkNode() const{
  return dynamic_cast<const FlatOCPInternal*>(get())!=0;
}

Variable& FlatOCP::variable(const std::string& name){
  return (*this)->variable(name);
}

void FlatOCP::addVariable(const std::string& name, const Variable& var){
  (*this)->addVariable(name,var);
}

void FlatOCP::makeAlgebraic(const std::string& name){
  (*this)->makeAlgebraic(variable(name));
}



} // namespace OptimalControl
} // namespace CasADi
