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
    mat_ = Matrix<double>(size1(),size2(),0);
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
  mat_ = Matrix<double>(size1(),size2(),col,rowind);
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

Matrix<double>& FunctionIO::mat(int dir){
  if(dir<0)
    return matA_.at(-1-dir);
  else if(dir==0)
    return mat_;
  else
    return matF_.at(-1+dir);
}

const Matrix<double>& FunctionIO::mat(int dir) const{
  if(dir<0)
    return matA_.at(-1-dir);
  else if(dir==0)
    return mat_;
  else
    return matF_.at(-1+dir);
}

Matrix<double>& FunctionIO::matF(int dir){
  return matF_.at(dir);
}

const Matrix<double>& FunctionIO::matF(int dir) const{
  return matF_.at(dir);
}

Matrix<double>& FunctionIO::matA(int dir){
  return matA_.at(dir);
}

const Matrix<double>& FunctionIO::matA(int dir) const{
  return matA_.at(dir);
}

vector<double>& FunctionIO::dataF(int dir){
  return matF(dir);
}

const vector<double>& FunctionIO::dataF(int dir) const{
  return matF(dir);
}

vector<double>& FunctionIO::dataA(int dir){
  return matA(dir);
}

const vector<double>& FunctionIO::dataA(int dir) const{
  return matA(dir);
}

int FunctionIO::numel() const{
  return size1()*size2();
}

int FunctionIO::size1() const{
  return mat_.size1();
}

int FunctionIO::size2() const{
  return mat_.size2();
}

void FunctionIO::getSparsityCRS(vector<int>& rowind, vector<int> &col) const{
  rowind = mat_.rowind();
  col = mat_.col();
}

void FunctionIO::getSparsity(vector<int>& row, vector<int> &col) const{
  col = mat_.col();
  row = mat_.sparsity().getRow();
}

int FunctionIO::size() const{
  return col().size();
}

int FunctionIO::sizeU() const{
  return mat_.sparsity().sizeU();
}

int FunctionIO::sizeL() const{
  return mat_.sparsity().sizeL();
}

vector<double>& FunctionIO::data(int dir){
  return mat(dir);
}

const vector<double>& FunctionIO::data(int dir) const{
  return mat(dir);
}

void FunctionIO::set(double val, int dir, Sparsity sp){
  mat(dir).set(val,sp);
}
    
void FunctionIO::get(double& val, int dir, Sparsity sp) const{
  mat(dir).get(val,sp);
}

void FunctionIO::set(const std::vector<double>& val, int dir, Sparsity sp){
  mat(dir).set(val,sp);
}

void FunctionIO::get(std::vector<double>& val, int dir, Sparsity sp) const{
  mat(dir).get(val,sp);
}

void FunctionIO::set(const double* val, int dir, Sparsity sp){
  mat(dir).set(val,sp);
}

void FunctionIO::get(double* val, int dir, Sparsity sp) const{
  mat(dir).get(val,sp);
}

const std::vector<int>& FunctionIO::rowind() const{
  return mat_.rowind();
}

std::vector<int>& FunctionIO::rowind(){
  return mat_.rowind();
}

const std::vector<int>& FunctionIO::col() const{
  return mat_.col();
}

std::vector<int>& FunctionIO::col(){
  return mat_.col();
}

int FunctionIO::rowind(int i) const{
  return mat_.rowind(i);
}

int FunctionIO::col(int el) const{
  return mat_.col(el);
}


} // namespace CasADi

