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
#include <cassert>
#include <typeinfo> 

using namespace std;

namespace CasADi{

MXNode::MXNode(){
  maxord_ = 0;
  nfdir_ = 1;
  nadir_ = 1;
}

MXNode::~MXNode(){
}

void MXNode::print(ostream &stream) const{
  stream << "<empty matrix expression>";
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

void MXNode::setOutput(const vector<double>& x){
  assert(x.size() == output().size());
  copy(x.begin(),x.end(), output().begin());
}

void MXNode::getOutput(vector<double>& x) const{
  assert(x.size() == output().size());
  copy(output().begin(),output().end(),x.begin());
}

void MXNode::setFwdSeed(const vector<double>& x, int dir){
  assert(x.size() == fwdSens(dir).size());
  copy(x.begin(),x.end(), fwdSens(dir).begin());
}

void MXNode::getFwdSens(vector<double>& x, int dir) const{
  assert(x.size() == fwdSens(dir).size());
  copy(fwdSens(dir).begin(),fwdSens(dir).end(),x.begin());
}

void MXNode::setAdjSeed(const vector<double>& x, int dir){
  assert(x.size() == adjSeed(dir).size());
  copy(x.begin(),x.end(), adjSeed(dir).begin());
}

void MXNode::getAdjSens(vector<double>& x, int dir) const{
  assert(x.size() == adjSeed(dir).size());
  copy(adjSeed(dir).begin(),adjSeed(dir).end(),x.begin());
}

MX& MXNode::dep(int ind){
  return dep_.at(ind);
}

const MX& MXNode::dep(int ind) const{
  return dep_.at(ind);
}
  
int MXNode::ndep() const{
  return dep_.size();
}

void MXNode::init(){
  forward_sensitivities_.resize(nfdir_);
  for(int dir=0; dir<nfdir_; ++dir)
    forward_sensitivities_[dir] = output_;
    
  adjoint_seeds_.resize(nadir_);
  for(int dir=0; dir<nadir_; ++dir)
    adjoint_seeds_[dir] = output_;
}

const std::vector<double>& MXNode::input(int ind) const{
  return dep_.at(ind)->output();
}

std::vector<double>& MXNode::input(int ind){
  return dep_.at(ind)->output();
}

const std::vector<double>& MXNode::output() const{
  return output_;
}

std::vector<double>& MXNode::output(){
  return output_;
}
  
const std::vector<double>& MXNode::fwdSeed(int ind, int dir) const{
  return dep_.at(ind)->fwdSens(dir);
}

std::vector<double>& MXNode::fwdSeed(int ind, int dir){
  return dep_.at(ind)->fwdSens(dir);
}

const std::vector<double>& MXNode::adjSeed(int dir) const{
  return adjoint_seeds_.at(dir);
}

std::vector<double>& MXNode::adjSeed(int dir){
  return adjoint_seeds_.at(dir);
}
  
const std::vector<double>& MXNode::fwdSens(int dir) const{
  return forward_sensitivities_.at(dir);
}

std::vector<double>& MXNode::fwdSens(int dir){
  return forward_sensitivities_.at(dir);
}

const std::vector<double>& MXNode::adjSens(int ind, int dir) const{
  return dep_.at(ind)->adjSeed(dir);
}

std::vector<double>& MXNode::adjSens(int ind, int dir){
  return dep_.at(ind)->adjSeed(dir);
}

void MXNode::setSize(int nrow, int ncol){
  output_ = Matrix<double>(nrow,ncol,0);
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

void MXNode::setDependencies(const std::vector<MX>& dep){
  dep_ = dep;
}

int MXNode::size1() const{
  return output_.size1();
}

int MXNode::size2() const{
  return output_.size2();
}



} // namespace CasADi
