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
#include "mx_tools.hpp"
#include "../fx/sx_function.hpp"
#include "evaluation_mx.hpp"
#include "symbolic_mx.hpp"
#include "constant_mx.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "norm.hpp"
#include "../casadi_math.hpp"

using namespace std;
namespace CasADi{

  MX::~MX(){
  }

  MX::MX(){
  }

  MX::MX(double x){
    assignNode(ConstantMX::create(sp_dense(1,1),x));
  }

  MX::MX(const Matrix<double>& x){
    assignNode(ConstantMX::create(x));
  }

  MX::MX(const vector<double>& x){
    assignNode(ConstantMX::create(x));
  }

  MX::MX(const string& name, int n, int m){
    assignNode(new SymbolicMX(name,n,m));
  }

  MX::MX(const std::string& name,const std::pair<int,int> &nm) {
    assignNode(new SymbolicMX(name,nm.first,nm.second));
  }

  MX::MX(const string& name, const CRSSparsity & sp){
    assignNode(new SymbolicMX(name,sp));
  }

  MX::MX(int nrow, int ncol){
    assignNode(new Constant<CompiletimeConst<0> >(CRSSparsity(nrow,ncol)));
  }

  MX::MX(const CRSSparsity& sp, const MX& val){
    // Make sure that val is dense and scalar
    casadi_assert(val.scalar());
  
    // Dense matrix if val dense
    if(val.dense()){
      *this = val->getGetNonzeros(sp,vector<int>(sp.size(),0));
     } else {
      // Empty matrix
      *this = sparse(sp.size1(),sp.size2());
    }
  }

  MX MX::create(MXNode* node){
    MX ret;
    ret.assignNode(node);
    return ret;
  }

  bool MX::__nonzero__() const {
    if (isNull()) {casadi_error("Cannot determine truth value of null MX.");}
    return (operator->())->__nonzero__();
  }

  const MX MX::sub(int i, const std::vector<int>& j) const{
    return sub(vector<int>(1,i),j);
  }

  const MX MX::sub(const std::vector<int>& i, int j) const{
    return sub(i,vector<int>(1,j));
  }

  const MX MX::sub(const vector<int>& ii, const vector<int>& jj) const{
    // Nonzero mapping from submatrix to full
    vector<int> mapping;
  
    // Get the sparsity pattern
    CRSSparsity sp = sparsity().sub(ii,jj,mapping);
 
    // Create return MX
    return (*this)->getGetNonzeros(sp,mapping);
  }

  const MX MX::sub(const Matrix<int>& k, int dummy) const {
    int s = size();
    const std::vector<int> & d = k.data();
    for (int i=0;i< k.size();i++) {
      if (d[i]>=s) casadi_error("MX::sub: a non-zero element at position " << d[i] << " was requested, but MX is only " << dimString());
    }
    return (*this)->getGetNonzeros(k.sparsity(),k.data());
  }

  const MX MX::sub(int i, int j) const{
    int ind = sparsity().getNZ(i,j);
    if (ind>=0) {
      return (*this)->getGetNonzeros(CRSSparsity::scalarSparsity,vector<int>(1,ind));
    } else {
      return (*this)->getGetNonzeros(CRSSparsity::scalarSparsitySparse,vector<int>(0));
    }
  }

  const MX MX::sub(const std::vector<int>& ii, const Matrix<int>& k) const{
    std::vector< int > cols = range(size2());
    std::vector< MX > temp;

    for (int i=0;i<ii.size();++i) {
      MX m(k.sparsity(),MX(0));
      for (int j=0;j<m.size();++j) {
        m[j] = sub(ii.at(i),k.at(j));
      }
      temp.push_back(m);
    }
    MX ret = vertcat(temp);
    simplify(ret);
    return ret;
  }

  const MX MX::sub(const Matrix<int>& k, const std::vector<int>& jj) const{
    std::vector< int > rows = range(size1());
    std::vector< MX > temp;

    for (int j=0;j<jj.size();++j) {
      MX m(k.sparsity(),MX(0));
      for (int i=0;i<m.size();++i) {
        m[i] = sub(k.at(i),jj.at(j));
      }
      temp.push_back(m);
    }
    MX ret = horzcat(temp);
    simplify(ret);
    return ret;
  }

  const MX MX::sub(const Matrix<int>& i, const Matrix<int>& j) const {
    casadi_assert_message(i.sparsity()==j.sparsity(),"sub(Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " and " << j.dimString() << ".");

    MX ret(i.sparsity(),MX(0));
    for (int k=0;k<i.size();++k) {
      ret[k] = sub(i.at(k),j.at(k));
    }
    simplify(ret);
    return ret;
  }

  const MX MX::sub(const CRSSparsity& sp, int dummy) const {
    casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"sub(CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
    vector<unsigned char> mappingc; // Mapping that will be filled by patternunion
  
    sparsity().patternCombine(sp, false, true, mappingc);
    vector<int> nz(sp.size(),-1);

    int k_this = 0;     // Non-zero of this matrix
    int k_sp = 0;       // Non-zero of resulting matrix
    for (vector<unsigned char>::const_iterator i=mappingc.begin(); i!=mappingc.end(); ++i){
      // In this matrix
      if(*i & 1){
        if(*i & 4){
          k_this++;
        } else {
          nz[k_sp++] = k_this++; // In both this matrix and in resulting matrix 
        }
      } else if(*i &2){
        k_sp++;
      }
    }

    MX ret = (*this)->getGetNonzeros(sp,nz);
    return ret;
  }

  void MX::setSub(const MX& m, int i, int j){
    setSub(m,vector<int>(1,i),vector<int>(1,j));
  }

  void MX::setSub(const MX& m, int i, const std::vector<int>& j){
    setSub(m,vector<int>(1,i),j);
  }

  void MX::setSub(const MX& m, const std::vector<int>& i, int j){
    setSub(m,i,vector<int>(1,j));
  }

  void MX::setSub(const MX& m, const Matrix<int>& k){
    // Allow m to be a 1x1
    if (m.dense() && m.scalar()) {
      if (k.numel()>1) {
        setSub(MX(k.sparsity(),m),k);
        return;
      }
    }
  
    casadi_assert_message(k.sparsity()==m.sparsity(),"Sparsity mismatch." << "lhs is " << k.dimString() << ", while rhs is " << m.dimString());
  
    casadi_error("MX::setSub not implemented yet");
  }

  void MX::setSub(const MX& m, const Slice& i, const Slice& j){
    setSub(m,i.getAll(size1()),j.getAll(size2()));
  }


  void MX::setSub(const MX& m, const vector<int>& ii, const vector<int>& jj){
    // Allow m to be a 1x1
    if (m.dense() && m.scalar()) {
      if (ii.size()>1 || jj.size()>1) {
        setSub(MX(ii.size(),jj.size(),m),ii,jj);
        return;
      }
    }
  
    casadi_assert_message(ii.size()==m.size1(),"Dimension mismatch." << "lhs is " << ii.size() << " x " << jj.size() << ", while rhs is " << m.dimString());
    casadi_assert_message(jj.size()==m.size2(),"Dimension mismatch." << "lhs is " << ii.size() << " x " << jj.size() << ", while rhs is " << m.dimString());

    if(dense() && m.dense()){
      // Dense mode
      int ld = size2(), ld_el = m.size2(); // leading dimensions
      vector<int> kk1, kk2;
      for(int i=0; i<ii.size(); ++i) {
        for(int j=0; j<jj.size(); ++j) {
          kk1.push_back(ii[i]*ld + jj[j]);
          kk2.push_back(i*ld_el+j);
        }
      }
      (*this)[kk1]=m[kk2];
    } else {
      // Sparse mode

      // Remove submatrix to be replaced
      erase(ii,jj);

      // Extend m to the same dimension as this
      MX el_ext = m;
      el_ext.enlarge(size1(),size2(),ii,jj);

      // Unite the sparsity patterns
      *this = unite(*this,el_ext);
    }
  }

  void MX::setSub(const MX& m, const Matrix<int>& i, const std::vector<int>& jj) {
    // If m is scalar
    if(m.scalar() && (jj.size() > 1 || i.size() > 1)){
      setSub(repmat(MX(i.sparsity(),m),1,jj.size()),i,jj);
      return;
    }

    if (!inBounds(jj,size2())) {
      casadi_error("setSub[.,i,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    //CRSSparsity result_sparsity = repmat(i,1,jj.size()).sparsity();
    CRSSparsity result_sparsity = horzcat(std::vector< Matrix<int> >(jj.size(),i)).sparsity();
  
    casadi_assert_message(result_sparsity == m.sparsity(),"setSub(.,Imatrix" << i.dimString() << ",Ivector(length=" << jj.size() << "),Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix,1," << jj.size() << ") = " << result_sparsity.dimString()  << " must match the sparsity of MX = "  << m.dimString() << ".");

    std::vector<int> slice_i = range(i.size1());
  
    for(int k=0; k<jj.size(); ++k) {
      MX el_k = m(slice_i,range(k*i.size2(),(k+1)*i.size2()));
      for (int j=0;j<i.size();++j) {
        (*this)(i.at(j),jj[k])=el_k[j];
      }
    }
  
  }



  void MX::setSub(const MX& m, const std::vector<int>& ii, const Matrix<int>& j) {
    // If m is scalar
    if(m.scalar() && (ii.size() > 1 || j.size() > 1)){
      setSub(repmat(MX(j.sparsity(),m),ii.size(),1),ii,j);
      return;
    }

    if (!inBounds(ii,size1())) {
      casadi_error("setSub[.,ii,j] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    //CRSSparsity result_sparsity = repmat(j,ii.size(),1).sparsity();
    CRSSparsity result_sparsity = vertcat(std::vector< Matrix<int> >(ii.size(),j)).sparsity();
  
    casadi_assert_message(result_sparsity == m.sparsity(),"setSub(Ivector(length=" << ii.size() << "),Imatrix" << j.dimString() << ",MX)::Dimension mismatch. The sparsity of repmat(Imatrix," << ii.size() << ",1) = " << result_sparsity.dimString() << " must match the sparsity of Matrix<T> = " << m.dimString() << ".");
  
    std::vector<int> slice_j = range(j.size2());
  
    for(int k=0; k<ii.size(); ++k) {
      MX el_k = m(range(k*j.size1(),(k+1)*j.size1()),slice_j);
      for (int i=0;i<j.size();++i) {
        (*this)(ii[k],j.at(i))=el_k[i];
      }
    }
  
  }


  void MX::setSub(const MX& m, const Matrix<int>& i, const Matrix<int>& j) {
    casadi_assert_message(i.sparsity()==j.sparsity(),"setSub(Imatrix m, Imatrix i, Imatrix j): sparsities must match. Got " << i.dimString() << " for i and " << j.dimString() << " for j.");

    // If m is scalar
    if(m.scalar() && i.numel() > 1){
      setSub(MX(i.sparsity(),m),i,j);
      return;
    }
  
    casadi_assert_message(m.sparsity()==i.sparsity(),"setSub(MX m, Imatrix i, Imatrix j): sparsities must match. Got " << m.dimString() << " for m and " << j.dimString() << " for i and j.");
  
    for(int k=0; k<i.size(); ++k) {
      (*this)(i.at(k),j.at(k)) = m[k]; 
    }
  }

  void MX::setSub(const MX& m, const CRSSparsity& sp, int dummy) {
    casadi_assert_message(size1()==sp.size1() && size2()==sp.size2(),"setSub(.,CRSSparsity sp): shape mismatch. This matrix has shape " << size1() << " x " << size2() << ", but supplied sparsity index has shape " << sp.size1() << " x " << sp.size2() << "." );
    casadi_error("Not implemented yet");
  }

  MX MX::getNZ(int k) const{
    if (k<0) k+=size();
    casadi_assert_message(k<size(),"MX::getNZ: requested at(" <<  k << "), but that is out of bounds:  " << dimString() << ".");
    return getNZ(vector<int>(1,k));
  }

  MX MX::getNZ(const vector<int>& k) const{
    CRSSparsity sp(k.size(),1,true);
  
    for (int i=0;i<k.size();i++) {
      casadi_assert_message(k[i] < size(),"Mapping::assign: index vector reaches " << k[i] << ", while dependant is only of size " << size());
    }
  
    MX ret = (*this)->getGetNonzeros(sp,k);
    return ret;
  }

  MX MX::getNZ(const Matrix<int>& k) const{
    CRSSparsity sp(k.size(),1,true);
    MX ret = (*this)->getGetNonzeros(sp,k.data());
    return ret;
  }

  void MX::setNZ(int k, const MX& el){
    if (k<0) k+=size();
    casadi_assert_message(k<size(),"MX::setNZ: requested at(" <<  k << "), but that is out of bounds:  " << dimString() << ".");
    setNZ(vector<int>(1,k),el);
  }

  void MX::setNZ(const vector<int>& k, const MX& el){
    casadi_assert_message(k.size()==el.size() || el.size()==1,
                          "MX::setNZ: length of non-zero indices (" << k.size() << ") " <<
                          "must match size of rhs (" << el.size() << ")."
                          );
  
    // Call recursively if points both objects point to the same node
    if(this==&el){
      MX el2 = el;
      setNZ(k,el2);
      return;
    }
  
    // Assert correctness
    for(int i=0; i<k.size(); ++i){
      casadi_assert_message(k[i] < size(), "Mapping::assign: index vector reaches " << k[i] << ", while dependant is only of size " << size());
    }
  
    // Quick return if no assignments to be made
    if(k.empty()) return;

    // Temporary
    MX x;

    // Project scalars
    if(k.size()!=el.size() && el.scalar() && el.dense()){      
      MX new_el = el->getGetNonzeros(sp_dense(k.size(),1),vector<int>(k.size(),0));
      x = new_el->getSetNonzeros(*this,k);
    } else {
      // Create a nonzero assignment node
      x = el->getSetNonzeros(*this,k);
    }
    *this = x;
    simplify(*this);
  }

  void MX::setNZ(const Matrix<int>& kk, const MX& m){
    if (m.size()==1 && m.numel()==1) {
      setNZ(kk.data(),m);
      return;
    }
    casadi_assert_message(kk.sparsity()==m.sparsity(),"Matrix<T>::setNZ: sparsity of IMatrix index " << kk.dimString() << " " << std::endl << "must match sparsity of rhs " << m.dimString() << ".");
    setNZ(kk.data(),m);
  }

  const MX MX::at(int k) const {
    return getNZ(k); 
  }

  /// Access a non-zero element
  NonZeros<MX,int> MX::at(int k) {
    return NonZeros<MX,int>(*this,k);
  }

  MX MX::binary(int op, const MX &x, const MX &y){
    return x->getBinarySwitch(op,y);
  }

  MX MX::unary(int op, const MX &x){
    return x->getUnary(Operation(op));
  }

  MXNode* MX::operator->(){
    return static_cast<MXNode*>(SharedObject::operator->());
  }

  const MXNode* MX::operator->() const{
    return static_cast<const MXNode*>(SharedObject::operator->());
  }

  MX MX::repmat(const MX& x, const std::pair<int, int> &nm){
    return repmat(x,nm.first,nm.second);
  }

  MX MX::repmat(const MX& x, int nrow, int ncol){
    if(x.scalar()){
      return MX(nrow,ncol,x);
    } else {
      casadi_assert_message(0,"not implemented");
      return MX();
    }
  }

  MX MX::sparse(int nrow, int ncol){
    return MX(nrow,ncol);
  }

  MX MX::sparse(const std::pair<int, int> &nm){
    return sparse(nm.first,nm.second);
  }

  MX MX::zeros(int nrow, int ncol){
    return zeros(sp_dense(nrow,ncol));
  }

  MX MX::zeros(const std::pair<int, int> &nm){
    return zeros(nm.first,nm.second);
  }

  MX MX::zeros(const CRSSparsity& sp){
    return create(ConstantMX::create(sp,0));
  }

  MX MX::ones(const CRSSparsity& sp){
    return create(ConstantMX::create(sp,1));
  }

  MX MX::ones(int nrow, int ncol){
    return ones(sp_dense(nrow,ncol));
  }

  MX MX::ones(const std::pair<int, int> &nm){
    return ones(nm.first,nm.second);
  }

  MX MX::inf(int nrow, int ncol){
    return inf(sp_dense(nrow,ncol));
  }

  MX MX::inf(const std::pair<int, int> &nm){
    return inf(nm.first,nm.second);
  }

  MX MX::inf(const CRSSparsity& sp){
    return create(ConstantMX::create(sp,numeric_limits<double>::infinity()));
  }

  MX MX::nan(int nrow, int ncol){
    return nan(sp_dense(nrow,ncol));
  }

  MX MX::nan(const std::pair<int, int> &nm){
    return nan(nm.first,nm.second);
  }

  MX MX::nan(const CRSSparsity& sp){
    return create(ConstantMX::create(sp,numeric_limits<double>::quiet_NaN()));
  }

  MX MX::eye(int n){
    Matrix<double> I(CRSSparsity::createDiagonal(n),1);
    return MX(I);
  }


  MX MX::operator-() const{

    if((*this)->getOp()==OP_NEG){
      return (*this)->dep(0);
    } else {
      return (*this)->getUnary(OP_NEG);
    }
  }

  MX::MX(const MX& x) : SharedObject(x){
  }

  const CRSSparsity& MX::sparsity() const{
    return (*this)->sparsity();
  }

  CRSSparsity& MX::sparsityRef(){
    // Since we can potentially change the behavior of the MX node, we must make a deep copy if there are other references
    makeUnique();
  
    // Return the reference, again, deep copy if multiple references
    (*this)->sparsity_.makeUnique();
    return (*this)->sparsity_;
  }

  void MX::erase(const vector<int>& ii, const vector<int>& jj){
    // Get sparsity of the new matrix
    CRSSparsity sp = sparsity();
  
    // Erase from sparsity pattern
    vector<int> mapping = sp.erase(ii,jj);
  
    // Create new matrix
    if(mapping.size()!=size()){
      MX ret = (*this)->getGetNonzeros(sp,mapping);
      *this = ret;
    }
  }

  void MX::enlarge(int nrow, int ncol, const vector<int>& ii, const vector<int>& jj){
    CRSSparsity sp = sparsity();
    sp.enlarge(nrow,ncol,ii,jj);
  
    MX ret = (*this)->getGetNonzeros(sp,range(size()));
    *this = ret;
  }

  MX::MX(int nrow, int ncol, const MX& val){
    // Make sure that val is scalar
    casadi_assert(val.scalar());
    casadi_assert(val.dense());
  
    CRSSparsity sp(nrow,ncol,true);
    *this = val->getGetNonzeros(sp,vector<int>(sp.size(),0));
  }

  MX MX::mul_full(const MX& y) const{
    const MX& x = *this;
    return x->getMultiplication(y);
  }

  MX MX::mul(const MX& y) const{
    return mul_smart(y);
  }

  MX MX::inner_prod(const MX& y) const{
    return (*this)->getInnerProd(y);
  }

  MX MX::outer_prod(const MX& y) const{
    return mul(trans(y));
  }

  MX MX::__pow__(const MX& n) const{
    if(n->getOp()==OP_CONST){
      return MX::binary(OP_CONSTPOW,*this,n);
    } else {
      return MX::binary(OP_POW,*this,n);
    }
  }

  MX MX::constpow(const MX& b) const{    
    return binary(OP_CONSTPOW,*this,b);
  }

  MX MX::fmin(const MX& b) const{       
    return binary(OP_FMIN,*this,b);
  }

  MX MX::fmax(const MX& b) const{       
    return binary(OP_FMAX,*this,b);
  }

  MX MX::arctan2(const MX& b) const{       
    return binary(OP_ATAN2,*this,b);
  }

  MX MX::printme(const MX& b) const{ 
    return binary(OP_PRINTME,*this,b);
  }

  MX MX::exp() const{ 
    return (*this)->getUnary(OP_EXP);
  }

  MX MX::log() const{ 
    return (*this)->getUnary(OP_LOG);
  }

  MX MX::log10() const{ 
    return log()*(1/std::log(10.));
  }

  MX MX::sqrt() const{ 
    return (*this)->getUnary(OP_SQRT);
  }

  MX MX::sin() const{ 
    return (*this)->getUnary(OP_SIN);
  }

  MX MX::cos() const{ 
    return (*this)->getUnary(OP_COS);
  }

  MX MX::tan() const{ 
    return (*this)->getUnary(OP_TAN);
  }

  MX MX::arcsin() const{ 
    return (*this)->getUnary(OP_ASIN);
  }

  MX MX::arccos() const{ 
    return (*this)->getUnary(OP_ACOS);
  }

  MX MX::arctan() const{ 
    return (*this)->getUnary(OP_ATAN);
  }

  MX MX::sinh() const{ 
    return (*this)->getUnary(OP_SINH);
  }

  MX MX::cosh() const{ 
    return (*this)->getUnary(OP_COSH);
  }

  MX MX::tanh() const{ 
    return (*this)->getUnary(OP_TANH);
  }

  MX MX::arcsinh() const{ 
    return (*this)->getUnary(OP_ASINH);
  }

  MX MX::arccosh() const{ 
    return (*this)->getUnary(OP_ACOSH);
  }

  MX MX::arctanh() const{ 
    return (*this)->getUnary(OP_ATANH);
  }

  MX MX::floor() const{ 
    return (*this)->getUnary(OP_FLOOR);
  }

  MX MX::ceil() const{ 
    return (*this)->getUnary(OP_CEIL);
  }

  MX MX::fabs() const{ 
    return (*this)->getUnary(OP_FABS);
  }

  MX MX::sign() const{ 
    return (*this)->getUnary(OP_SIGN);
  }

  MX MX::erfinv() const{ 
    return (*this)->getUnary(OP_ERFINV);
  }

  MX MX::erf() const{ 
    return (*this)->getUnary(OP_ERF);
  }

  MX MX::logic_not() const{ 
    return (*this)->getUnary(OP_NOT);
  }

  void MX::lift(const MX& x_guess){ 
    casadi_assert(sparsity()==x_guess.sparsity());
    *this = (*this)->getBinary(OP_LIFT,x_guess,false,false);
  }

  MX MX::__add__(const MX& y) const{
    return MX::binary(OP_ADD,*this,y);
  }

  MX MX::__sub__(const MX& y) const{
    return MX::binary(OP_SUB,*this,y);
  }

  MX MX::__mul__(const MX& y) const{
    return MX::binary(OP_MUL,*this,y);
  }

  MX MX::__div__(const MX& y) const{
    return MX::binary(OP_DIV,*this,y);
  }

  MX MX::__lt__(const MX& y) const{
    return MX::binary(OP_LT,*this,y);
  }
  
  MX MX::__le__(const MX& y) const{
    return MX::binary(OP_LE,*this,y);
  }
  
  MX MX::__eq__(const MX& y) const{
    return MX::binary(OP_EQ,*this,y);
  }
  
  MX MX::__ne__(const MX& y) const{
    return MX::binary(OP_NE,*this,y);
  }
  
  MX MX::logic_and(const MX& y) const{
    return MX::binary(OP_AND,*this,y);
  }
  
  MX MX::logic_or(const MX& y) const{
    return MX::binary(OP_OR,*this,y);
  }
  
  MX MX::if_else_zero(const MX& y) const{
    return MX::binary(OP_IF_ELSE_ZERO,*this,y);
  }
  
  MX MX::__constpow__(const MX& b) const { return (*this).constpow(b);}
  MX MX::__mrdivide__(const MX& b) const { if (b.scalar()) return *this/b; throw CasadiException("mrdivide: Not implemented");}
  MX MX::__mpower__(const MX& b) const   { return pow(*this,b); throw CasadiException("mpower: Not implemented");}

  void MX::append(const MX& y){
    *this = vertcat(*this,y);
  }

  long MX::max_num_calls_in_print_ = 10000;

  void MX::setMaxNumCallsInPrint(long num){
    max_num_calls_in_print_ = num;
  }

  long MX::getMaxNumCallsInPrint(){
    return max_num_calls_in_print_;
  }

  MX MX::getDep(int ch) const { return !isNull() ? (*this)->dep(ch) : MX(); }

  int MX::getNdeps() const { return !isNull() ? (*this)->ndep() : 0; }
  
  std::string MX::getName() const { return !isNull() ? (*this)->getName() : "null"; }

  bool         MX::isSymbolic () const { return !isNull() ? (*this)->getOp()==OP_PARAMETER : false; }
  bool         MX::isConstant () const { return !isNull() ? (*this)->getOp()==OP_CONST : false; }
  bool         MX::isEvaluation () const { return !isNull() ? (*this)->getOp()==OP_CALL : false; }
  bool         MX::isEvaluationOutput () const { return !isNull() ? (*this)->isOutputNode() : false; }

  int         MX::getEvaluationOutput () const { return !isNull() ? (*this)->getFunctionOutput() : -1; }

    
  bool         MX::isOperation (int op) const { return !isNull() ? (*this)->getOp()==op : false; }
  bool         MX::isMultiplication () const { return !isNull() ? (*this)->getOp()==OP_MATMUL : false; }

  bool         MX::isNorm () const { return !isNull() ? dynamic_cast<const Norm*>(get())!=0 : false; }

  FX MX::getFunction () {  return (*this)->getFunction(); }
         
  double MX::getValue() const{
    return (*this)->getValue();
  }

  Matrix<double> MX::getMatrixValue() const{
    return (*this)->getMatrixValue();
  }
           
  bool MX::isBinary() const { return isNull() ? false : (*this)->isBinaryOp();}

  bool MX::isUnary() const { return isNull() ? false : (*this)->isUnaryOp();}
         
  int MX::getOp() const {
    return (*this)->getOp();
  }

  bool MX::isEqual(const MX& y, int depth) const{
    return isEqual(static_cast<const MXNode*>(y.get()),depth);
  }

  bool MX::isEqual(const MXNode* y, int depth) const{
    if(get()==y)
      return true;
    else if(depth>0)
      return (*this)->isEqual(y,depth);
    else
      return false;
  }

  bool MX::isCommutative() const {
    if (isUnary()) return true;
    casadi_assert_message(isBinary() || isUnary(),"MX::isCommutative: must be binary or unary operation");
    return operation_checker<CommChecker>(getOp());
  }

  long MX::__hash__() const {
    if (isNull()) return 0;
    // FIXME: This is bad coding, the pointer might on certain architectures be larger than the long,
    // giving an error of type "error: cast from 'const SharedObjectInternal*' to 'int' loses precision"
    // The solution is to convert to intptr_t, but this is not yet C++ standard
    return long(get());
  }

  Matrix<int> MX::mapping() const {
    return (*this)->mapping();
  }

  int MX::getTemp() const{
    return (*this)->temp;
  }
    
  void MX::setTemp(int t){
    (*this)->temp = t;
  }

  int MX::getNumOutputs() const{
    return (*this)->getNumOutputs();
  }
  
  MX MX::getOutput(int oind) const{
    return (*this)->getOutput(oind);
  }

  MX MX::setSparse(const CRSSparsity& sp) const{
    if(sp==sparsity()){
      return *this;
    } else {
      return (*this)->getSetSparse(sp);
    }
  }

 MX MX::makeDense(const MX& val) const{
    if(dense()){
      return *this;
    } else {
      MX ret(size1(),size2(),val);
      ret(sparsity()) = *this;
      return ret;
    }
  }

  int MX::eq_depth_ = 1;

  void MX::setEqualityCheckingDepth(int eq_depth){
    eq_depth_ = eq_depth;
  }

  int MX::getEqualityCheckingDepth(){
    return eq_depth_;
  }

          
} // namespace CasADi
