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

#include "fx_internal.hpp"
#include "../mx/evaluation.hpp"
#include <typeinfo> 
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"

using namespace std;

namespace CasADi{
  
FX::FX(){
}

FX::~FX(){
}

const FXInternal* FX::operator->() const{
  return (const FXInternal*)OptionsFunctionality::operator->();
}

FXInternal* FX::operator->(){
  return (FXInternal*)OptionsFunctionality::operator->();
}

MX FX::operator()(const MX &x, int ind) const{
  vector<MX> dep(1);
  dep[0] = x;
  MX ret;
  ret.assignNode(new Evaluation(*this,dep,ind));
  return ret;
}

MX FX::operator()(const vector<MX> &x, int ind) const{
  MX ret;
  ret.assignNode(new Evaluation(*this,x,ind));
  return ret;
}

void FX::evaluate(int fsens_order, int asens_order){
  (*this)->assertInit();
  (*this)->evaluate(fsens_order,asens_order);
}

void FX::solve(){
  (*this)->evaluate(0,0);
}

int FX::getNumInputs() const{
  return (*this)->input_.size();
}

int FX::getNumOutputs() const{
  return (*this)->output_.size();  
}

void FX::setNumInputs(int num_in){
  return (*this)->input_.resize(num_in);
}

void FX::setNumOutputs(int num_out){
  return (*this)->output_.resize(num_out);  
}


FX FX::jacobian(int iind, int oind){
  return (*this)->jacobian(iind,oind);  
}

FX FX::hessian(int iind, int oind){
  return (*this)->hessian(iind,oind);  
}

bool FX::isInit() const {
  return (*this)->is_init_;
}

void FX::assertInit() const{
  if(!(*this)->is_init_)
    throw CasadiException("FX::assertInit: function has not been initialized");
}

bool FX::checkNode() const{
  return dynamic_cast<const FXInternal*>(get());
}

FunctionIO& FX::input(int i){
  return (*this)->input(i);
}

const FunctionIO& FX::input(int i) const{
  return (*this)->input(i);
}
  
FunctionIO& FX::output(int i){
  return (*this)->output(i);
}

const FunctionIO& FX::output(int i) const{
  return (*this)->output(i);
}

const vector<double>& FX::getInputData(int ind) const {
  return input(ind).data(); 
}

const vector<double>& FX::getOutputData(int ind) const {
  return output(ind).data(); 
}

const vector<double>& FX::getFwdSeedData(int ind, int dir) const {return input(ind).dataF(dir); }
const vector<double>& FX::getFwdSensData(int ind, int dir) const {return output(ind).dataF(dir); }
const vector<double>& FX::getAdjSeedData(int ind, int dir) const {return output(ind).dataA(dir); }
const vector<double>& FX::getAdjSensData(int ind, int dir) const {return input(ind).dataA(dir); }
vector<double>& FX::getInputData(int ind) {return input(ind).data(); }
vector<double>& FX::getOutputData(int ind) {return output(ind).data(); }
vector<double>& FX::getFwdSeedData(int ind, int dir) {return input(ind).dataF(dir); }
vector<double>& FX::getFwdSensData(int ind, int dir) {return output(ind).dataF(dir); }
vector<double>& FX::getAdjSeedData(int ind, int dir) {return output(ind).dataA(dir); }
vector<double>& FX::getAdjSensData(int ind, int dir) {return input(ind).dataA(dir); }






} // namespace CasADi

