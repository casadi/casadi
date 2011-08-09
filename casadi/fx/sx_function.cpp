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

#include "sx_function_internal.hpp"
#include <cassert>
#include <limits>
#include <stack>
#include <deque>
#include <fstream>
#include <sstream>
#include "../stl_vector_tools.hpp"
#include "../sx/sx_node.hpp"
#include "../mx/evaluation.hpp"
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"

namespace CasADi{

using namespace std;


SXFunction::SXFunction(){
}

SXFunction::SXFunction(const SXMatrix& arg, const SXMatrix& res){
  vector<SXMatrix> argv(1,arg);
  vector<SXMatrix> resv(1,res);
  assignNode(new SXFunctionInternal(argv,resv));
}

SXFunction::SXFunction(const vector< vector<SX> >& arg, const vector< vector<SX> >& res){
  vector<SXMatrix> argv(arg.begin(),arg.end());
  vector<SXMatrix> resv(res.begin(),res.end());
  assignNode(new SXFunctionInternal(argv,resv));
}

SXFunction::SXFunction(const vector< SXMatrix>& arg, const vector<SXMatrix>& res){
  assignNode(new SXFunctionInternal(arg,res));
}

SXFunction::SXFunction(const vector< vector<SX> >& arg, const SXMatrix& res){
  vector<SXMatrix> argv(arg.begin(),arg.end());
  vector<SXMatrix> resv(1,res);
  assignNode(new SXFunctionInternal(argv,resv));
}

SXFunction::SXFunction(const vector< SXMatrix>& arg, const SXMatrix& res){
  vector<SXMatrix> resv(1,res);
  assignNode(new SXFunctionInternal(arg,resv));
}

const SXFunctionInternal* SXFunction::operator->() const{
  return (const SXFunctionInternal*)FX::operator->();
}

SXFunctionInternal* SXFunction::operator->(){
  return (SXFunctionInternal*)FX::operator->();
}

SXFunction SXFunction::jacobian(int iind, int oind){
  return shared_cast<SXFunction>(FX::jacobian(iind,oind));  
}

SXFunction SXFunction::hessian(int iind, int oind){
  return shared_cast<SXFunction>(FX::hessian(iind,oind));
}

bool SXFunction::checkNode() const{
  return dynamic_cast<const SXFunctionInternal*>(get())!=0;
}

SXMatrix SXFunction::jac(int iind, int oind){
  if(input(iind).empty() || output(oind).empty()) return Matrix<SX>(); // quick return
  vector<pair<int,int> > jblocks(1,pair<int,int>(oind,iind));
  if(!isInit()) init();
  return jac(jblocks).front();
}

SXMatrix SXFunction::grad(int iind, int oind){
  return trans(jac(iind,oind));
}

SXMatrix SXFunction::hess(int iind, int oind){
  if(output(oind).numel() != 1)
    throw CasadiException("SXFunctionInternal::hess: function must be scalar");
  
  // Reverse mode to calculate gradient
  if((*this)->verbose()){
    cout << "SXFunction::hess: calculating gradient " << endl;
  }
  Matrix<SX> g = grad(iind,oind);
  if((*this)->verbose()){
    cout << "SXFunction::hess: calculating gradient done " << endl;
  }
  
  // Create function
  SXFunction gfcn(inputSX(iind),g);

  // Make gradient verbose
  if((*this)->verbose()){
    gfcn.setOption("verbose",true);
  }
  
  // Initialize
  gfcn.init();
  
  // Calculate jacobian of gradient
  if((*this)->verbose()){
    cout << "SXFunction::hess: calculating Jacobian " << endl;
  }
  SXMatrix ret = gfcn.jac();
  if((*this)->verbose()){
    cout << "SXFunction::hess: calculating Jacobian done" << endl;
  }
  
  // Return jacobian of the gradient
  return ret;
}

const SXMatrix& SXFunction::inputSX(int ind) const{
  return (*this)->inputv_[ind];
}

const SXMatrix& SXFunction::outputSX(int ind) const{
  return (*this)->outputv_[ind];
}

void SXFunction::generateCode(const string& filename){
  (*this)->generateCode(filename);
}

const vector<SXAlgEl>& SXFunction::algorithm() const{
  return (*this)->algorithm_;
}

int SXFunction::countNodes() const{
  casadi_assert(isInit());
  return algorithm().size();
}

void SXFunction::clearSymbolic(){
  (*this)->clearSymbolic();
}

vector<Matrix<SX> > SXFunction::jac(const vector<pair<int,int> >& jblocks){
  return (*this)->jac(jblocks);
}

SXFunction::SXFunction(const MXFunction& f){
  MXFunction f2 = f;
  SXFunction t = f2.expand();
  assignNode(t.get());
}



} // namespace CasADi

