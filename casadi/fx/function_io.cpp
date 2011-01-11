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

#include "function_io.hpp"
#include <sstream>
#include "../stl_vector_tools.hpp"
#include "../casadi_exception.hpp"

using namespace std;

namespace CasADi{

FunctionIO::FunctionIO(){
  dense_ = true;
}

void FunctionIO::init(){
  
  // Make dense if necessary
  if(dense_){
    get() = Matrix<double>(get().size1(),get().size2(),0);
  }

  // Non-zeros
  for(int i=0; i<matF_.size(); ++i)
    matF_[i] = mat_;
  for(int i=0; i<matA_.size(); ++i)
    matA_[i] = mat_;
}

void FunctionIO::setSize(int nrow, int ncol){
  mat_.resize(nrow,ncol);
}

void FunctionIO::setSparsityCRS(const vector<int>& rowind, const vector<int> &col){
  mat_ = Matrix<double>(get().size1(),get().size2(),col,rowind);
  dense_ = false;
}

void FunctionIO::setNumFwdDir(int nfdir){
  matF_.resize(nfdir);
}

void FunctionIO::setNumAdjDir(int nadir){
  matA_.resize(nadir);
}

int FunctionIO::numFwdDir() const{
  return matF_.size();
}

int FunctionIO::numAdjDir() const{
  return matA_.size();
}

Matrix<double>& FunctionIO::get(){
  return mat_;
}
    
const Matrix<double>& FunctionIO::get() const{
  return mat_;
}
    
Matrix<double>& FunctionIO::getFwd(int dir){
  return matF_.at(dir);
}
    
const Matrix<double>& FunctionIO::getFwd(int dir) const{
  return matF_.at(dir);
}
    
Matrix<double>& FunctionIO::getAdj(int dir){
  return matA_.at(dir);
}
    
const Matrix<double>& FunctionIO::getAdj(int dir) const{
  return matA_.at(dir);
}


} // namespace CasADi

