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
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"

namespace CasADi{

using namespace std;


SXFunction::SXFunction(){
}

SXFunction::SXFunction(const SXMatrix& arg, const SXMatrix& res){
  assignNode(new SXFunctionInternal(vector<SXMatrix>(1,arg),
                                    vector<SXMatrix>(1,res)));
}

SXFunction::SXFunction(const vector< vector<SX> >& arg, const vector< vector<SX> >& res){
  assignNode(new SXFunctionInternal(vector<SXMatrix>(arg.begin(),arg.end()),
                                    vector<SXMatrix>(res.begin(),res.end())));
}

SXFunction::SXFunction(const vector< SXMatrix>& arg, const vector<SXMatrix>& res){
  assignNode(new SXFunctionInternal(arg,res));
}

SXFunction::SXFunction(const vector< SXMatrix>& arg, const IOSchemeVector< SXMatrix >& res){
  assignNode(new SXFunctionInternal(arg,res));
  setOutputScheme(res.io_scheme());
}

SXFunction::SXFunction(const IOSchemeVector< SXMatrix >& arg, const vector< SXMatrix>& res){
  assignNode(new SXFunctionInternal(arg,res));
  setInputScheme(arg.io_scheme());
}

SXFunction::SXFunction(const IOSchemeVector< SXMatrix >& arg, const IOSchemeVector< SXMatrix >& res){
  assignNode(new SXFunctionInternal(arg,res));
  setInputScheme(arg.io_scheme());
  setOutputScheme(res.io_scheme());
}

SXFunction::SXFunction(const vector< vector<SX> >& arg, const SXMatrix& res){
  assignNode(new SXFunctionInternal(vector<SXMatrix>(arg.begin(),arg.end()),
                                    vector<SXMatrix>(1,res)));
}

SXFunction::SXFunction(const vector< SXMatrix>& arg, const SXMatrix& res){
  assignNode(new SXFunctionInternal(arg,
                                    vector<SXMatrix>(1,res)));
}

SXFunction::SXFunction(const SXMatrix& arg, const std::vector< std::vector<SX> >& res){
  assignNode(new SXFunctionInternal(vector<SXMatrix>(1,arg),
                                    vector<SXMatrix>(res.begin(),res.end())));
  
}

SXFunction::SXFunction(const SXMatrix& arg, const std::vector< SXMatrix>& res){
  assignNode(new SXFunctionInternal(vector<SXMatrix>(1,arg),
                                    res));
}


const SXFunctionInternal* SXFunction::operator->() const{
  return static_cast<const SXFunctionInternal*>(FX::operator->());
}

SXFunctionInternal* SXFunction::operator->(){
  return static_cast<SXFunctionInternal*>(FX::operator->());
}

bool SXFunction::checkNode() const{
  return dynamic_cast<const SXFunctionInternal*>(get())!=0;
}

SXMatrix SXFunction::jac(int iind, int oind, bool compact, bool symmetric){
  return (*this)->jac(iind,oind,compact,symmetric);
}

SXMatrix SXFunction::grad(int iind, int oind){
  return (*this)->grad(iind,oind);
}

SXMatrix SXFunction::hess(int iind, int oind){
  return (*this)->hess(iind,oind);
}

const SXMatrix& SXFunction::inputExpr(int ind) const{
  return (*this)->inputv_.at(ind);
}

const SXMatrix& SXFunction::outputExpr(int ind) const{
  return (*this)->outputv_.at(ind);
}
  
const std::vector<SXMatrix>& SXFunction::inputExpr() const{
  return (*this)->inputv_;
}
  
const std::vector<SXMatrix> & SXFunction::outputExpr() const{
  return (*this)->outputv_;
}

const vector<ScalarAtomic>& SXFunction::algorithm() const{
  return (*this)->algorithm_;
}

int SXFunction::countNodes() const{
  assertInit();
  return algorithm().size() - getNumScalarOutputs();
}

void SXFunction::clearSymbolic(){
  (*this)->clearSymbolic();
}

SXFunction::SXFunction(const MXFunction& f){
  MXFunction f2 = f;
  SXFunction t = f2.expand();
  assignNode(t.get());
}

SXFunction::SXFunction(const FX& f){
  const SXFunctionInternal* temp = dynamic_cast<const SXFunctionInternal*>(f.get());
  if (temp) {
    assignNode(temp->clone());
  } else {
    MXFunction f2(f);
    SXFunction t = f2.expand();
    assignNode(t.get());
  }
}

SXFunction SXFunction::operator[](int k) const {

  // Delegate to FX
  MXFunction temp = shared_cast<MXFunction>(shared_cast<FX>(*this)[k]);
  
  casadi_assert(!temp.isNull());
  
  // Expand to SXFunction
  SXFunction ret = temp.expand(inputExpr());

  ret.init();

  return ret;
}

std::vector<SX> SXFunction::getFree() const{
  return (*this)->free_vars_;
}

int SXFunction::getWorkSize() const{
  return (*this)->work_.size();
}

} // namespace CasADi

