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
#include "../fx/mx_function.hpp"
#include <typeinfo> 
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"
#include "parallelizer.hpp"

using namespace std;

namespace CasADi{
  
FX::FX(){
}

FX::~FX(){
}

FX FX::create(FXInternal* node){
  FX ret;
  ret.assignNode(node);
  return ret;
}

const FXInternal* FX::operator->() const{
  return (const FXInternal*)OptionsFunctionality::operator->();
}

FXInternal* FX::operator->(){
  return (FXInternal*)OptionsFunctionality::operator->();
}

vector<MX> FX::call(const MX &x){
  vector<MX> xvec(1,x);
  return call(xvec);
}

vector<MX> FX::call(const vector<MX> &x){
  assertInit();
  
  MX ev;
  ev.assignNode(new Evaluation(*this,x));
  vector<MX> ret(getNumOutputs());
  for(int i=0; i<ret.size(); ++i){
    if(output(i).numel()>0)
      ret[i].assignNode(new OutputNode(ev,i));
  }
  return ret;
}

vector<vector<MX> > FX::call(const vector<vector<MX> > &x, const Dictionary& paropt){
  assertInit();
  
  // Make sure not empty
  casadi_assert_message(x.size()>1,"FX: call(vector<vector<MX> >): argument must be of length > 1. You supplied length " << x.size() << ".");
  
  // Return object
  vector<vector<MX> > ret(x.size());
  

  
  // Check if we are bypassing the parallelizer
  Dictionary::const_iterator ii=paropt.find("parallelization");
  if(ii!=paropt.end() && ii->second=="expand"){
    for(int i=0; i<x.size(); ++i){
      ret[i] = call(x[i]);
    }
    return ret;
  }
    
//  if(paropt
    
//     getOption("parallelization")=="serial")

  
  // Create parallelizer object and initialize it
  Parallelizer p(vector<FX>(x.size(),*this));
  p.setOption(paropt);
  p.init();
  
  // Concatenate the arguments
  vector<MX> p_in;
  p_in.reserve(x.size() * getNumInputs());
  for(int i=0; i<x.size(); ++i){
    p_in.insert(p_in.end(),x[i].begin(),x[i].end());
  }
  
  // Call the parallelizer
  vector<MX> p_out = p.call(p_in);
  casadi_assert(p_out.size() == x.size() * getNumOutputs());

  // Collect the outputs
  vector<MX>::const_iterator it=p_out.begin();
  for(int i=0; i<x.size(); ++i){
    ret[i].insert(ret[i].end(),it,it+getNumOutputs());
    it += getNumOutputs();
  }
  return ret;
}

void FX::evaluate(int nfdir, int nadir){
  assertInit();
  casadi_assert(nfdir<=(*this)->nfdir_);
  casadi_assert(nadir<=(*this)->nadir_);
  (*this)->evaluate_switch(nfdir,nadir);
}

void FX::evaluate_old(int fsens_order, int asens_order){
  evaluate(fsens_order * (*this)->nfdir_, 
           asens_order * (*this)->nadir_);
}

void FX::solve(){
  evaluate(0,0);
}

int FX::getNumInputs() const{
  return (*this)->getNumInputs();
}

int FX::getNumOutputs() const{
  return (*this)->getNumOutputs();  
}

void FX::setNumInputs(int num_in){
  return (*this)->setNumInputs(num_in);
}

void FX::setNumOutputs(int num_out){
  return (*this)->setNumOutputs(num_out);  
}

FX FX::jacobian(int iind, int oind){
  assertInit();
  casadi_assert(iind>=0 && iind<getNumInputs());
  casadi_assert(oind>=0 && oind<getNumOutputs());

  vector<pair<int,int> > jblocks;
  jblocks.push_back(pair<int,int>(oind,iind));
  
  if((*this)->store_jacobians_){
    // Get a reference to the place the Jacobian is or will be saved
    FX& J = (*this)->jacs_[iind][oind];
    
    // Generate a Jacobian if necessary
    if(J.isNull())
      J = (*this)->jacobian_switch(jblocks);
    
    // Return a reference to the stored Jacobian
    return J;
    
  } else {
    return (*this)->jacobian_switch(jblocks);
  }
}

FX FX::jacobian(const std::vector<std::pair<int,int> >& jblocks){
  assertInit();
  return (*this)->jacobian_switch(jblocks);
}

FX FX::hessian(int iind, int oind){
  assertInit();
  return (*this)->hessian(iind,oind);  
}

bool FX::checkNode() const{
  return dynamic_cast<const FXInternal*>(get())!=0;
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

void FX::addMonitor(const string& mon){
  (*this)->monitors_.insert(mon);
}

void FX::removeMonitor(const string& mon){
  (*this)->monitors_.erase(mon);
}

const Dictionary & FX::getStats() const{
  return (*this)->getStats();
}

GenericType FX::getStat(const string& name) const{
  return (*this)->getStat(name);
}

CRSSparsity& FX::jacSparsity(int iind, int oind){
  return (*this)->jacSparsity(iind,oind);
}

void FX::setJacSparsity(const CRSSparsity& sp, int iind, int oind){
  (*this)->setJacSparsity(sp,iind,oind);
}

std::vector<MX> FX::symbolicInput() const{
  return (*this)->symbolicInput();
}

std::vector<SXMatrix> FX::symbolicInputSX() const{
  return (*this)->symbolicInputSX();
}

FX FX::operator[](int k) const {

  // Argument checking
  if (k<0) k+=getNumOutputs();
  casadi_assert_message(k<getNumOutputs(),"FX[int k]:: Attempt to select the k'th output with k=" << k << ", but should be smaller or equal to number of outputs (" << getNumOutputs() << ").");
  
  // Get the inputs in MX form
  std::vector< MX > in = symbolicInput();
  
  // Clone such that we can use const.
  FX clone = *this;
  
  // Get the outputs
  std::vector< MX > result = clone.call(in);
  
  // Construct an MXFunction with only the k'th output
  MXFunction ret(in,result[k]);
  
  // Initialize it
  ret.init();
  
  // And return it, will automatically cast to FX
  return ret;
}

#if 0

vector<DMatrix> FX::jac(const vector<DMatrix> &x, int iind){
  casadi_assert(0);  
}

vector<SXMatrix> FX::jac(const vector<SXMatrix> &x, int iind){
  casadi_assert(0);
}

vector<MX> FX::jac(const vector<MX> &x, int iind){
  casadi_assert(0);
}

vector<DMatrix> FX::jac(const vector<DMatrix> &x, const vector<DMatrix> &v){
  casadi_assert(0);
}

vector<SXMatrix> FX::jac(const vector<SXMatrix> &x, const vector<SXMatrix> &v){
  casadi_assert(0);
}

vector<MX> FX::jac(const vector<MX> &x, const vector<MX> &v){
  casadi_assert(0);  
}

vector<DMatrix> FX::grad(const vector<DMatrix> &x, int oind){
    casadi_assert(0);
}

vector<SXMatrix> FX::grad(const vector<SXMatrix> &x, int oind){
  casadi_assert(0);  
}

vector<MX> FX::grad(const vector<MX> &x, int oind){
  casadi_assert(0);  
}

vector<DMatrix> FX::grad(const vector<DMatrix> &x, const vector<DMatrix> &v){
  casadi_assert(0);  
}

vector<SXMatrix> FX::grad(const vector<SXMatrix> &x, const vector<SXMatrix> &v){
  casadi_assert(0);
}

vector<MX> FX::grad(const vector<MX> &x, const vector<MX> &v){
  casadi_assert(0);
}

#endif

} // namespace CasADi

