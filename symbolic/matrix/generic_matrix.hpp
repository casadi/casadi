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

#ifndef GENERIC_MATRIX_HPP
#define GENERIC_MATRIX_HPP

#include "slice.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "crs_sparsity.hpp"
#include "sparsity_tools.hpp"
#include "../casadi_math.hpp"

namespace CasADi{

  /** Sparsity format for getting and setting inputs and outputs */
  enum Sparsity{SPARSE,SPARSESYM,DENSE,DENSESYM};

  /** \brief Matrix base class
  This is a common base class for MX and Matrix<>, introducing a uniform syntax and implementing
  common functionality using the curiously recurring template pattern (CRTP) idiom.\n
  
  The class is designed with the idea that "everything is a matrix", that is, also scalars and vectors.\n
  This philosophy makes it easy to use and to interface in particularily with Python and Matlab/Octave.\n
  
  The syntax tries to stay as close as possible to the ublas syntax  when it comes to vector/matrix operations.\n

  Index starts with 0.\n
  Index flatten happens as follows: (i,j) -> k = j+i*size2()\n
  Vectors are considered to be column vectors.\n
  
  The storage format is a (modified) compressed row storage (CRS) format. This way, a vector element can always be accessed in constant time.\n
  
  The sparsity can be accessed with CRSSparsity& sparsity()\n
  
  \author Joel Andersson 
  \date 2012    
*/
template<typename MatType>
class GenericMatrix{
  public:
    
    /** \brief Get the number of (structural) non-zero elements */
    int size() const;

    /** \brief Get the number of non-zeros in the lower triangular half */
    int sizeL() const;

    /** \brief Get get the number of non-zeros in the upper triangular half */
    int sizeU() const;

    /** \brief Get get the number of non-zeros on the diagonal */
    int sizeD() const;
    
    /** \brief Get the number of elements */
    int numel() const;

    /** \brief Get the first dimension (i.e. n for a n-by-m matrix) */
    int size1() const;
    
    /** \brief Get the first dimension (i.e. m for a n-by-m matrix) */
    int size2() const;

    /** \brief Get the number if non-zeros for a given sparsity pattern */
    int size(Sparsity sp) const;
    
    /** \brief Get string representation of dimensions.
    The representation is (nrow x ncol = numel | size)
    */
    std::string dimString() const;
    
    #ifndef SWIG  
    /** \brief  Get the shape */
    std::pair<int,int> shape() const;
    #endif

    /** \brief  Check if the matrix expression is empty, i.e. one of its dimensions is 0 */
    bool empty() const;
    
    /** \brief  Check if the matrix expression is null, i.e. its dimensios are 0-by-0 */
    bool null() const;
    
    /** \brief  Check if the matrix expression is dense */
    bool dense() const;
    
    /** \brief  Check if the matrix expression is scalar */
    bool scalar(bool scalar_and_dense=false) const;

    /** \brief  Check if the matrix expression is square */
    bool square() const;

    /** \brief Get the sparsity pattern */
    const CRSSparsity& sparsity() const;

    /** \brief Access the sparsity, make a copy if there are multiple references to it */
    CRSSparsity& sparsityRef();

    #ifndef SWIG
    /** \brief  Get vector nonzero or slice of nonzeros */
    template<typename K>
    const MatType operator[](const K& k) const{ return static_cast<const MatType*>(this)->getNZ(k); }

    /** \brief  Access vector nonzero or slice of nonzeros */
    template<typename K>
    NonZeros<MatType,K> operator[](const K& k){ return NonZeros<MatType,K>(static_cast<MatType&>(*this),k); }

    /** \brief  Get vector element or slice */
    template<typename I>
    const MatType operator()(const I& i) const{ return static_cast<const MatType*>(this)->sub(i,0);}

    /** \brief  Get Sparsity slice */
    const MatType operator()(const CRSSparsity& sp) const{ return static_cast<const MatType*>(this)->sub(sp); }
    
    /** \brief  Get Matrix element or slice */
    template<typename I, typename J>
    const MatType operator()(const I& i, const J& j) const{ return static_cast<const MatType*>(this)->sub(i,j); }

    /** \brief  Access vector element or slice */
    template<typename I>
    SubMatrix<MatType,I,int> operator()(const I& i){ return SubMatrix<MatType,I,int>(static_cast<MatType&>(*this),i,0); }

    /** \brief  Access Sparsity slice */
    SubMatrix<MatType,CRSSparsity,int> operator()(const CRSSparsity& sp){ return SubMatrix<MatType,CRSSparsity,int>(static_cast<MatType&>(*this),sp,0); }
      
    /** \brief  Access Matrix element or slice */
    template<typename I, typename J>
    SubMatrix<MatType,I,J> operator()(const I& i, const J& j){ return SubMatrix<MatType,I,J>(static_cast<MatType&>(*this),i,j); }
    #endif // SWIG

    /** \brief Create an n-by-m matrix with symbolic variables */
    static MatType sym(const std::string& name, int n=1, int m=1);

    /** \brief Create a vector of length p with with matrices with symbolic variables of given sparsity */
    static std::vector<MatType > sym(const std::string& name, const CRSSparsity& sp, int p);

    /** \brief Create a vector of length p with n-by-m matrices with symbolic variables */
    static std::vector<MatType > sym(const std::string& name, int n, int m, int p);

    /** \brief Create an matrix with symbolic variables, given a sparsity pattern */
    static MatType sym(const std::string& name, const CRSSparsity& sp);
    
    /** \brief Matrix-matrix multiplication.
    * Attempts to identify quick returns on matrix-level and 
    * delegates to MatType::mul_full if no such quick returns are found.
    */
    MatType mul_smart(const MatType& y, const CRSSparsity& sp_z) const;
    
};

#ifndef SWIG
// Implementations

template<typename MatType>
const CRSSparsity& GenericMatrix<MatType>::sparsity() const{
  return static_cast<const MatType*>(this)->sparsity();
}

template<typename MatType>
CRSSparsity& GenericMatrix<MatType>::sparsityRef(){
  return static_cast<MatType*>(this)->sparsityRef();
}

template<typename MatType>
int GenericMatrix<MatType>::size() const{
  return sparsity().size();
}

template<typename MatType>
int GenericMatrix<MatType>::sizeU() const{
  return sparsity().sizeU();
}

template<typename MatType>
int GenericMatrix<MatType>::sizeL() const{
  return sparsity().sizeL();
}

template<typename MatType>
int GenericMatrix<MatType>::sizeD() const{
  return sparsity().sizeD();
}

template<typename MatType>
int GenericMatrix<MatType>::numel() const{
  return sparsity().numel();
}

template<typename MatType>
int GenericMatrix<MatType>::size1() const{
  return sparsity().size1();
}

template<typename MatType>
int GenericMatrix<MatType>::size2() const{
  return sparsity().size2();
}

template<typename MatType>
std::pair<int,int> GenericMatrix<MatType>::shape() const{
  return sparsity().shape();
}

template<typename MatType>
std::string GenericMatrix<MatType>::dimString() const {
  return sparsity().dimString();
}

template<typename MatType>
bool GenericMatrix<MatType>::empty() const{
  return numel()==0;
}

template<typename MatType>
bool GenericMatrix<MatType>::null() const{
  return size1()==0 && size2()==0;
}

template<typename MatType>
bool GenericMatrix<MatType>::dense() const{
  return numel()==size();
}

template<typename MatType>
bool GenericMatrix<MatType>::scalar(bool scalar_and_dense) const{
  return sparsity().scalar(scalar_and_dense);
}

template<typename MatType>
bool GenericMatrix<MatType>::square() const{
  return sparsity().square();
}

template<typename MatType>
MatType GenericMatrix<MatType>::mul_smart(const MatType& y, const CRSSparsity &sp_z) const {
  const MatType& x = *static_cast<const MatType*>(this);
  
  if (!(x.scalar() || y.scalar())) {
    casadi_assert_message(size2()==y.size1(),"Matrix product with incompatible dimensions. Lhs is " << dimString() << " and rhs is " << y.dimString() << ".");
  }
  
  // Check if we can simplify the product
  if(isIdentity(x)){
    return y;
  } else if(isIdentity(y)){
    return x;
  } else if(isZero(x) || isZero(y)){
    // See if one of the arguments can be used as result
    if(x.size()==0 && y.size1()==y.size2()) {
      return x;
    } else if(y.size()==0 && x.size1()==x.size2()) {
      return y;
    } else {
      if (x.size()==0 || y.size()==0 || y.empty() || x.empty()) {
        return MatType::sparse(x.size1(),y.size2());
      } else {
        return MatType::zeros(x.size1(),y.size2());
      }
    }
  } else if(x.scalar() || y.scalar()){
    return x*y;
  } else {
    return x.mul_full(y,sp_z);
  }
}

template<typename MatType>
int GenericMatrix<MatType>::size(Sparsity sp) const{
  if(sp==SPARSE){
    return size();
  } else if(sp==SPARSESYM){
    return sizeL();
  } else if(sp==DENSE){
    return numel();
  } else if(sp==DENSESYM){
    return (numel()+size1())/2;
  } else {
      throw CasadiException("Matrix<T>::size(Sparsity): unknown sparsity");
  }
}


#endif // SWIG

template<typename MatType>
MatType GenericMatrix<MatType>::sym(const std::string& name, int n, int m){ return sym(name,sp_dense(n,m));}

template<typename MatType>
std::vector<MatType> GenericMatrix<MatType>::sym(const std::string& name, const CRSSparsity& sp, int p){
  std::vector<MatType> ret(p);
  std::stringstream ss;
  for(int k=0; k<p; ++k){
    ss.str("");
    ss << name << k;
    ret[k] = sym(ss.str(),sp);
  }
  return ret;
}

template<typename MatType>
std::vector<MatType > GenericMatrix<MatType>::sym(const std::string& name, int n, int m, int p){ return sym(name,sp_dense(n,m),p);}

template<typename MatType>
MatType GenericMatrix<MatType>::sym(const std::string& name, const CRSSparsity& sp){ throw CasadiException("\"sym\" not defined for instantiation");}

} // namespace CasADi

#endif // GENERIC_MATRIX_HPP

