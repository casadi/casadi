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

#include "symbolic_mx_node.hpp"
// #include "../sx/sx_tools.hpp"

using namespace std;

namespace CasADi{

SymbolicMatrix::SymbolicMatrix(const std::string& name, int n, int m) : name_(name) {
  setSparsity(CRSSparsity(n,m,true));
/*  sv_ = symbolic(name,n,m);*/
}

SymbolicMatrix::SymbolicMatrix(const std::string& name, const CRSSparsity & sp) : name_(name) {
  setSparsity(sp);
/*  sv_ = symbolic(name,n,m);*/
}

SymbolicMatrix* SymbolicMatrix::clone() const{
  return new SymbolicMatrix(*this);
}

void SymbolicMatrix::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << name_;
}

void SymbolicMatrix::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
}

void SymbolicMatrix::evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output){
}

bool SymbolicMatrix::isSymbolic() const{
  return true;
}

const std::string& SymbolicMatrix::getName() const{
  return name_;
}

std::vector<MX> SymbolicMatrix::partial(const std::vector<MX>& x){
  return std::vector<MX>(1,MX::eye(sparsity().numel()));
}


} // namespace CasADi

