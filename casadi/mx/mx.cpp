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

#include "mx.hpp"
#include "mx_node.hpp"
#include "../sx/sx_matrix.hpp"
#include "scalar_matrix_op.hpp"
#include "matrix_scalar_op.hpp"
#include "matrix_matrix_op.hpp"
#include "unary_op.hpp"
#include "if_else_node.hpp"
#include "matrix_element.hpp"
#include "../fx/sx_function.hpp"
#include "evaluation.hpp"
#include "symbolic_mx_node.hpp"
#include "mx_constant.hpp"
#include "slice.hpp"

namespace CasADi{

MX::~MX(){
}

MX::MX(){
}

MX::MX(double x){
  assignNode(new MXConstant(&x,1,1));
}

MX::MX(const std::vector<double> &x){
  assignNode(new MXConstant(&x[0],x.size(),1,'R'));  
}

MX::MX(const std::vector<double> &x, int n, int m, char order){
  if(x.size() != n*m) throw CasadiException("MX(const std::vector<double> &, ...): dimension mismatch");
  assignNode(new MXConstant(&x[0],n,m,order));  
}

MX::MX(const std::string& name, int n, int m){
  assignNode(new SymbolicMatrix(name,n,m));
}

MX::Element MX::operator()(int i, int j){
  return Element(*this, j+i*size2());
}

MX::Element MX::operator[](int k){
 return Element(*this, k);
}

MX MX::getElement(int k) const{
  MX ret;
  ret.assignNode(new MatrixElement(*this,k/size2(), k%size2()));
  return ret;  
}

MX MX::slice(Slicer i, Slicer j) const{
  MX ret;
  ret.assignNode(new Slice(*this,i,j));
  return ret;
}

MX MX::getRow(int i) const {return slice(i,SlicerPrimitiveAll());}

MX MX::getColumn(int j) const {return slice(SlicerPrimitiveAll(),j);}

MX& MX::setElement(const MX& el, int k){
  throw CasadiException("MX::setElement: not implemented");
}

const MX MX::operator()(int i, int j) const{
  return getElement(j+i*size2());
}

const MX MX::operator[](int k) const{
  return getElement(k);
}

const MatrixSize& MX::size() const{
  return (*this)->sz;
}

int MX::size1() const{
  return (*this)->sz.nrow;
}

int MX::numel() const{
  return size1()*size2();
}

int MX::size2() const{
  return (*this)->sz.ncol;
}

MX operator+(const MX &x, const MX &y){
  return MX::binary(ADD_NODE,x,y);
}

MX operator-(const MX &x, const MX &y){
  return MX::binary(SUB_NODE,x,y);
}

MX MX::binary(int op, const MX &x, const MX &y){
  if(x.numel()==1)
    return scalar_matrix(op,x,y);
  else if(y.numel()==1)  
    return matrix_scalar(op,x,y);
  else
    return matrix_matrix(op,x,y);
}

MX MX::unary(int op, const MX &x){
  MX ret;
  ret.assignNode(new UnaryOp(OPERATION(op),x));
  return ret;
}

MX MX::scalar_matrix(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new ScalarMatrixOp(OPERATION(op),x,y));  
  return ret;
}

MX MX::matrix_scalar(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new MatrixScalarOp(OPERATION(op),x,y));  
  return ret;
}

MX MX::matrix_matrix(int op, const MX &x, const MX &y){
  MX ret;
  ret.assignNode(new MatrixMatrixOp(OPERATION(op),x,y)); 
  return ret;
}

MX operator*(const MX &x, const MX &y){
  return MX::binary(MUL_NODE,x,y);
}

MX operator/(const MX &x, const MX &y){
  return MX::binary(DIV_NODE,x,y);
}

MXNode* MX::operator->(){
  return (MXNode*)SharedObject::operator->();
}

const MXNode* MX::operator->() const{
  return (const MXNode*)SharedObject::operator->();
}

MXNode* MX::get(){
  return (MXNode*)SharedObject::get();
}

const MXNode* MX::get() const{
  return (const MXNode*)SharedObject::get();
}

MX& MX::operator+=(const MX &y){
  return *this = *this + y;
}

MX& MX::operator-=(const MX &y){
  return *this = *this - y;
}

MX& MX::operator*=(const MX &y){
  return *this = *this * y;
}

MX& MX::operator/=(const MX &y){
  return *this = *this / y;
}

MX if_else(const MX &cond, const MX &if_true, const MX &if_false){
  MX ret;
  ret.assignNode(new IfElseNode(cond,if_true,if_false));
  return ret;
}

MX::Element::Element(MX& mx_, int k_) : mx(mx_), k(k_){ 
}

MX::Element::operator MX() const{
  return mx.getElement(k);
}

MX& MX::Element::operator=(const MX &y){
   return mx.setElement(y,k);
}

MX& MX::Element::operator+=(const MX &y){
   return mx.setElement(mx.getElement(k)+y,k);
}

MX& MX::Element::operator-=(const MX &y){
   return mx.setElement(mx.getElement(k)-y,k);  
}

MX& MX::Element::operator*=(const MX &y){
     return mx.setElement(mx.getElement(k)*y,k);
}

MX& MX::Element::operator/=(const MX &y){
    return mx.setElement(mx.getElement(k)/y,k); 
}

bool MX::isEmpty() const{
  return numel()==0;
}

MX::MX(const MatrixSize& sz){
  vector<double> v(sz.nrow*sz.ncol,0);
  assignNode(new MXConstant(&v[0],sz.nrow,sz.ncol,'R'));
}

MX MX::zeros(int nrow, int ncol){
  vector<double> v(nrow*ncol,0);
  return MX(v,nrow,ncol);
}

MX MX::ones(int nrow, int ncol){
  vector<double> v(nrow*ncol,1);
  return MX(v,nrow,ncol);
}

MX MX::operator-() const{
  return unary(NEG_NODE,*this);
}

void MX::Element::print(std::ostream &stream) const{
  mx.getElement(k).print(stream);
}

} // namespace CasADi


using namespace CasADi;
namespace std{

MX exp(const MX& x){
  return MX::unary(EXP_NODE,x);
}

MX log(const MX& x){
  return MX::unary(LOG_NODE,x);
}

MX sqrt(const MX& x){
  return MX::unary(SQRT_NODE,x);
}

MX sin(const MX& x){
  return MX::unary(SIN_NODE,x);
}

MX cos(const MX& x){;  
  return MX::unary(COS_NODE,x);
}

MX tan(const MX& x){
  return MX::unary(TAN_NODE,x);
}

MX atan(const MX& x){
  return MX::unary(ATAN_NODE,x);
}

MX asin(const MX& x){
  return MX::unary(ASIN_NODE,x);
}

MX acos(const MX& x){
  return MX::unary(ACOS_NODE,x);
}

MX pow(const MX& x, const MX& n){
  return MX::binary(POW_NODE,x,n);
}

MX erf(const MX& x){
  return MX::unary(ERF_NODE,x);
}

MX floor(const MX& x){
  return MX::unary(FLOOR_NODE,x);
}

MX ceil(const MX& x){
  return MX::unary(CEIL_NODE,x);
}

MX fmin(const MX& x, const MX& y){
  return MX::binary(FMIN_NODE,x,y);
}

MX fmax(const MX& x, const MX& y){
  return MX::binary(FMAX_NODE,x,y);
}



} // namespace std
