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

#include "mx_node.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <typeinfo> 
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace CasADi{

MXNode::MXNode(){
  temp = 0;
}

MXNode::~MXNode(){
}

const string& MXNode::getName() const{
  throw CasadiException(string("MXNode::getName() not defined for class ") + typeid(*this).name());
}

bool MXNode::isSymbolic() const{
  return false;
}

bool MXNode::isConstant() const{
  return false;
}

const MX& MXNode::dep(int ind) const{
  return dep_.at(ind);
}
  
MX& MXNode::dep(int ind){
  return dep_.at(ind);
}
  
int MXNode::ndep() const{
  return dep_.size();
}

void MXNode::setSparsity(const CRSSparsity& sparsity){
  sparsity_ = sparsity;
}

void MXNode::setDependencies(const MX& dep){
  dep_.resize(1);
  dep_[0] = dep;
}
    
void MXNode::setDependencies(const MX& dep1, const MX& dep2){
  dep_.resize(2);
  dep_[0] = dep1;
  dep_[1] = dep2;
}
    
void MXNode::setDependencies(const MX& dep1, const MX& dep2, const MX& dep3){
  dep_.resize(3);
  dep_[0] = dep1;
  dep_[1] = dep2;
  dep_[2] = dep3;
}

int MXNode::addDependency(const MX& dep){
  dep_.push_back(dep);
  return dep_.size()-1;
}

void MXNode::addDependency(int depind, const std::vector<int>& nz_d, const std::vector<int>& nz){
  casadi_assert(0);
}

void MXNode::addDependency(const MX& d, const std::vector<int>& nz_d, const std::vector<int>& nz){
  casadi_assert(0);
}
    
void MXNode::addDependency(const MX& d, const std::vector<int>& nz_d){
  casadi_assert(0);
}

void MXNode::setDependencies(const std::vector<MX>& dep){
  dep_ = dep;
}

int MXNode::numel() const{
  return sparsity_.numel();
}

int MXNode::size() const{
  return sparsity_.size();
}

int MXNode::size1() const{
  return sparsity_.size1();
}

int MXNode::size2() const{
  return sparsity_.size2();
}

const CRSSparsity& MXNode::sparsity() const{
  return sparsity_;
}

const CRSSparsity& MXNode::sparsity(int oind){
  casadi_assert_message(oind==0, "Index out of bounds");
  return sparsity_;
}

void MXNode::print(std::ostream &stream) const{
  vector<string> args(ndep());
  for(int i=0; i<ndep(); ++i){
    stringstream ss;
    if (dep(i).isNull()) {
      args[i] = "MX()";
      continue;
    } 
    dep(i)->print(ss);
    args[i] = ss.str();
  }
  print(stream,args);
}

FX& MXNode::getFunction(){
  throw CasadiException(string("MXNode::getFunction() not defined for class ") + typeid(*this).name());
}

int MXNode::getFunctionOutput() const{
  throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") + typeid(*this).name());
}

int MXNode::getFunctionInput() const{
  throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") + typeid(*this).name());
}

void MXNode::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output){
  DMatrixPtrVV fwdSeed, fwdSens, adjSeed, adjSens;
  evaluate(input,output,fwdSeed, fwdSens, adjSeed, adjSens);
}

void MXNode::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output){
  SXMatrixPtrVV fwdSeed, fwdSens, adjSeed, adjSens;
  evaluateSX(input,output,fwdSeed, fwdSens, adjSeed, adjSens);
}

void MXNode::evaluateMX(const MXPtrV& input, MXPtrV& output){
  MXPtrVV fwdSeed, fwdSens, adjSeed, adjSens;
  evaluateMX(input,output,fwdSeed, fwdSens, adjSeed, adjSens,false);
}

void MXNode::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  SharedObjectNode::deepCopyMembers(already_copied);
  dep_ = deepcopy(dep_,already_copied);
}




} // namespace CasADi
