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

vector<MX> FX::call(const MX &x) const{
  vector<MX> xvec(1,x);
  return call(xvec);
}

vector<MX> FX::call(const vector<MX> &x) const{
  vector<MX> ret(getNumOutputs());
  for(int i=0; i<ret.size(); ++i){
    if(output(i).numel()>0)
      ret[i].assignNode(new Evaluation(*this,x,i));
  }
  return ret;
}

#ifndef USE_FUNCTORS
MX FX::operator()(const MX &x, int ind) const{
  vector<MX> dep(1,x);
  MX ret;
  ret.assignNode(new Evaluation(*this,dep,ind));
  return ret;
}

MX FX::operator()(const vector<MX> &x, int ind) const{
  MX ret;
  ret.assignNode(new Evaluation(*this,x,ind));
  return ret;
}
#endif // USE_FUNCTORS

void FX::evaluate(int fsens_order, int asens_order){
  casadi_assert(isInit());
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
  return (*this)->setNumInputs(num_in);
}

void FX::setNumOutputs(int num_out){
  return (*this)->setNumOutputs(num_out);  
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

bool FX::checkNode() const{
  return dynamic_cast<const FXInternal*>(get());
}

Matrix<double>& FX::input(int iind){
  return (*this)->input(iind);
}
    
const Matrix<double>& FX::input(int iind) const{
  return (*this)->input(iind);
}

Matrix<double>& FX::output(int oind){
  return (*this)->output(oind);
}
    
const Matrix<double>& FX::output(int oind) const{
  return (*this)->output(oind);
}

Matrix<double>& FX::fwdSeed(int iind, int dir){
  return (*this)->fwdSeed(iind,dir);
}
    
const Matrix<double>& FX::fwdSeed(int iind, int dir) const{
  return (*this)->fwdSeed(iind,dir);
}

Matrix<double>& FX::fwdSens(int oind, int dir){
  return (*this)->fwdSens(oind,dir);
}
    
const Matrix<double>& FX::fwdSens(int oind, int dir) const{
  return (*this)->fwdSens(oind,dir);
}

Matrix<double>& FX::adjSeed(int oind, int dir){
  return (*this)->adjSeed(oind,dir);
}
    
const Matrix<double>& FX::adjSeed(int oind, int dir) const{
  return (*this)->adjSeed(oind,dir);
}

Matrix<double>& FX::adjSens(int iind, int dir){
  return (*this)->adjSens(iind,dir);
}
    
const Matrix<double>& FX::adjSens(int iind, int dir) const{
  return (*this)->adjSens(iind,dir);
}

void FX::addMonitor(const std::string& mon){
  (*this)->monitors_.insert(mon);
}

void FX::removeMonitor(const std::string& mon){
  (*this)->monitors_.erase(mon);
}



} // namespace CasADi

