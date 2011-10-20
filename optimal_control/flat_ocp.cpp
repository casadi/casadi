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
#include "variable_tools.hpp"

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

SX FlatOCP::t() const{
  return (*this)->t_;
}

std::vector<Variable>& FlatOCP::x(){
  return (*this)->x_;
}

std::vector<Variable>& FlatOCP::xd(){
  return (*this)->xd_;
}

std::vector<Variable>& FlatOCP::xa(){
  return (*this)->xa_;
}

std::vector<Variable>& FlatOCP::q(){
  return (*this)->q_;
}

std::vector<Variable>& FlatOCP::y(){
  return (*this)->y_;
}

std::vector<Variable>& FlatOCP::p(){
  return (*this)->p_;
}

std::vector<Variable>& FlatOCP::u(){
  return (*this)->u_;
}

std::vector<SX>& FlatOCP::dae(){
  return (*this)->dae_;
}
    
std::vector<SX>& FlatOCP::ode(){
  return (*this)->ode_;
}
    
std::vector<SX>& FlatOCP::alg(){
  return (*this)->alg_;
}
    
std::vector<SX>& FlatOCP::quad(){
  return (*this)->quad_;
}
    
std::vector<SX>& FlatOCP::dep(){
  return (*this)->dep_;
}
    
std::vector<SX>& FlatOCP::initial(){
  return (*this)->initial_;
}

std::vector<SX>& FlatOCP::mterm(){
  return (*this)->mterm_;
}

std::vector<SX>& FlatOCP::lterm(){
  return (*this)->lterm_;
}

std::vector<SX>& FlatOCP::path(){
  return (*this)->path_;
}

std::vector<double>& FlatOCP::path_min(){
  return (*this)->path_min_;
}

std::vector<double>& FlatOCP::path_max(){
  return (*this)->path_max_;
}

double& FlatOCP::t0(){
  return (*this)->t0_;
}

bool& FlatOCP::t0_free(){
  return (*this)->t0_free_;
}

double& FlatOCP::tf(){
  return (*this)->tf_;
}

bool& FlatOCP::tf_free(){
  return (*this)->tf_free_;
}

void FlatOCP::eliminateDependent(bool eliminate_dependents_with_bounds){
  (*this)->eliminateDependent(eliminate_dependents_with_bounds);
}

void FlatOCP::eliminateInterdependencies(){
  (*this)->eliminateInterdependencies();
}


} // namespace OptimalControl
} // namespace CasADi
