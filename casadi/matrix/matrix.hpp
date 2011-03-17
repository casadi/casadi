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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include "../casadi_exception.hpp"
#include "../printable_object.hpp"
#include "../casadi_limits.hpp"
#include "../casadi_operators.hpp"
#include "element.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "crs_sparsity.hpp"

namespace CasADi{
  
  /** Sparsity format for getting and setting inputs and outputs */
  enum Sparsity{SPARSE,SPARSESYM,DENSE,DENSESYM};

  /// Dummy class denoting all rows/columns
  class AllRange{
  };
  static AllRange ALL;

  /** \brief General sparse matrix class
  General sparse matrix class that is designed with the idea that "everything is a matrix", that is, also scalars and vectors.\n
  This philosophy makes it easy to use and to interface in particularily with Matlab and Python.\n
  
  The syntax tries to stay as close as possible to the ublas syntax  when it comes to vector/matrix operations.\n

  Index starts with 0.\n
  Index flatten happens as follows: (i,j) -> k = j+i*size2()\n
  Vectors are considered to be column vectors.\n
  
  The storage format is a (modified) compressed row storage (CRS) format. This way, a vector element can always be accessed in constant time.\n
  
  Matrix<T> is polymorphic with a std::vector<T> that contain all non-identical-zero elements.\n
  The sparsity can be accessed with CRSSparsity& sparsity()\n
  
  \author Joel Andersson 
  \date 2010	
*/
template<class T>
class Matrix : public std::vector<T>, public PrintableObject{

  public:
    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    Matrix();
    
    /// Copy constructor (normal)
    Matrix(const Matrix<T>& m);
    
#ifndef SWIG
    /// Assignment (normal)
    Matrix<T>& operator=(const Matrix<T>& m);
#endif // SWIG
    
    /// empty n-by-m matrix constructor
    Matrix(int n, int m);
    
    /// dense n-by-m matrix filled with val constructor
    Matrix(int n, int m, const T& val);

    /// sparse n-by-m matrix filled with given sparsity
    Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind, const std::vector<T>& data=std::vector<T>());

    /// sparse matrix with a given sparsity
    explicit Matrix(const CRSSparsity& sparsity, const T& val=0);
    
    /// This constructor enables implicit type conversion from a numeric type
    Matrix(double val);

    /// Construct from a vector
    /**
    * Thanks to implicit conversion, you can pretend that Matrix(const SX& x); exists.
    * Note: above remark applies only to C++, not python or octave interfaces
    */
    Matrix(const std::vector<T>& x);
    
    /// Construct dense matrix from a vector with the elements in column major ordering
    Matrix(const std::vector<T>& x, int n, int m);

#ifndef SWIG
    /// This constructor enables implicit type conversion from a matrix element
    Matrix(const Element<Matrix<T>,T>& val);

    /// Expose iterators
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::reverse_iterator reverse_iterator;
    typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
    
    /// References
    typedef T& reference;
    typedef const T& const_reference;
    
    /// Get iterators to beginning and end
    iterator begin(){ return data().begin();}
    const_iterator begin() const{ return data().begin();}
    reverse_iterator rbegin(){ return data().rbegin();}
    const_reverse_iterator rbegin() const{ return data().rbegin();}
    iterator end(){ return data().end();}
    const_iterator end() const{ return data().end();}
    reverse_iterator rend(){ return data().rend();}
    const_reverse_iterator rend() const{ return data().rend();}
    
    /// Get references to beginning and end
    reference front(){ return data().front();}
    const_reference front() const { return data().front();}
    reference back(){ return data().back();}
    const_reference back() const { return data().back();}
    
    /** \brief  Create an expression from an stl vector  */
    template<typename A>
    Matrix(const std::vector<A>& x){
      sparsity_ = CRSSparsity(x.size(),1,true);
      data().resize(x.size());
      copy(x.begin(),x.end(),begin());
    }

    /** \brief  Create a non-vector expression from an stl vector */
    template<typename A>
    Matrix(const std::vector<A>& x,  int n, int m){
      if(x.size() != n*m) throw CasadiException("Matrix::Matrix(const std::vector<T>& x,  int n, int m): dimension mismatch");
      sparsity_ = CRSSparsity(n,m,true);
      data().resize(x.size());
      copy(x.begin(),x.end(),begin());
    }
    
    /** \brief  ublas vector */
#ifdef HAVE_UBLAS
    template<typename T, typename A>
    explicit Matrix<T>(const ublas::vector<A> &x){
      sparsity_ = CRSSparsity(x.size(),1,true);
      std::vector<T>::resize(x.size());
      copy(x.begin(),x.end(),begin());
    }

    template<typename T, typename A>
    explicit Matrix<T>(const ublas::matrix<A> &x){
      sparsity_ = CRSSparsity(x.size1(),x.size2(),true);
      std::vector<T>::resize(numel());
      copy(x.begin(),x.end(),begin());
      return ret;
    }
#endif // HAVE_UBLAS
#endif // SWIG

    /// get the number of non-zeros
    int size() const;

    /// get the number of non-zeros in the lower triangular half
    int sizeL() const;

    /// get the number of non-zeros in the upper triangular half
    int sizeU() const;

    /// get the number of elements
    int numel() const;

    /// get the first dimension
    int size1() const;       
    
    /// get the second dimension
    int size2() const;

    //@{
    /// Check type of matrix
    bool empty() const; // is the matrix empty
    bool scalar() const; // is the matrix scalar
    bool vector() const; // is the matrix a vector
    bool dense() const; // is the matrix dense
    //@}


#ifndef SWIG    
    /// get an element
    const T getElement(int i=0, int j=0) const;
    
    /// set an element
    void setElement(int i, int j, const T& el);

    /// get a reference to an element
    T& getElementRef(int i=0, int j=0);
  
    /// Get a submatrix
    Matrix<T> getSub(const std::vector<int>& ii, const std::vector<int>& jj) const;
    
    /// Set a submatrix
    void setSub(const std::vector<int>& ii, const std::vector<int>& jj, const Matrix<T>& m);
  
    /// Get nonzeros
    Matrix<T> getNZ(const std::vector<int>& kk) const;
    
    /// Set nonzeros
    void setNZ(const std::vector<int>& kk, const Matrix<T>& m);
    
    /// Access an element 
    Element<Matrix<T>,T> operator()(int i, int j=0){ return Element<Matrix<T>,T>(*this,i,j); }

    /// Get an element 
    const T operator()(int i, int j=0) const{ return getElement(i,j); }

    /// Access a submatrix
    SubMatrix<Matrix<T> > operator()(const std::vector<int>& ii, const std::vector<int>& jj){ return SubMatrix<Matrix<T> >(*this,ii,jj);}
    
    /// Get a submatrix
    const Matrix<T> operator()(const std::vector<int>& ii, const std::vector<int>& jj) const{ return getSub(ii,jj);}

    /// Access a row
    SubMatrix<Matrix<T> > operator()(int i, const std::vector<int>& jj){ return SubMatrix<Matrix<T> >(*this,std::vector<int>(1,i),jj);}
    
    /// Get a row
    const Matrix<T> operator()(int i, const std::vector<int>& jj) const{ return getSub(std::vector<int>(1,i),jj);}

    /// Access a column
    SubMatrix<Matrix<T> > operator()(const std::vector<int>& ii, int j){ return SubMatrix<Matrix<T> >(*this,ii,std::vector<int>(1,j));}
    
    /// Get a column
    const Matrix<T> operator()(const std::vector<int>& ii, int j) const{ return getSub(ii,std::vector<int>(1,j));}

    /// Access all rows
    template<class A>
    const Matrix<T> operator()(const AllRange& i, A j) const{ return operator()(range(size1()),j);}
    
    /// Get all rows
    template<class A>
    SubMatrix<Matrix<T> > operator()(const AllRange& i, A j){ return operator()(range(size1()),j);}
  
    /// Access all columns
    template<class A>
    const Matrix<T> operator()(A i, const AllRange& j) const{ return operator()(i,range(size2()));}

    /// Get all columns
    template<class A>
    SubMatrix<Matrix<T> > operator()(A i, const AllRange& j){ return operator()(i,range(size2()));}

    /// Get all rows and columns
    SubMatrix<Matrix<T> > operator()(const AllRange& i, const AllRange& j){ return operator()(range(size1()),range(size2()));}

    /// Get a non-zero element
    const T& at(int k) const{ return data().at(k); }

    /// Access a non-zero element
    T& at(int k){ return data().at(k); }
    
    /// Get a non-zero element
    const T& operator[](int k) const{ return data().operator[](k); }

    /// Access a non-zero element
    T& operator[](int k){ return data().operator[](k); }
    
    /// Get a number of non-zero elements
    const Matrix<T> operator[](const std::vector<int>& kk) const{ return getNZ(kk);}
    
    /// Access a number of non-zero elements
    NonZeros<Matrix<T> > operator[](const std::vector<int>& kk){ return NonZeros<Matrix<T> >(*this,kk);}

#endif // SWIG

    /// Python: get a non-zero entry
    const T getitem(int i) const;
    
    /// Python: get a matrix entry
    const T getitem(const std::vector<int> &I) const;
    
    /// Python: get a submatrix
    const Matrix<T> getitem(const std::vector< std::vector<int> > &II) const;
    
    /// Python: set a non-zero entry
    void setitem(int k, const T& el);
    
    /// Python: set a matrix entry
    void setitem(const std::vector<int> &I, const T&  el);

    /// Python: set a submatrix
    void setitem(const std::vector< std::vector<int> > &II, const Matrix<T>& m);
    
    /// Set all elements to zero
    void setZero();
    
    /// Set all elements to a value
    void setAll(const T& val);
    
    /** \brief  Make the matrix an dense n-by-m matrix */
    void makeDense(int n, int m, const T& val);

    /** \brief  Make the matrix an empty n-by-m matrix */
    void makeEmpty(int n, int m);

    Matrix<T> operator+() const;
    Matrix<T> operator-() const;

    /** \brief  Unary function */
#ifndef SWIG
    Matrix<T> unary(T (*fcn)(const T&)) const;
    Matrix<T> binary(T (*fcn)(const T&, const T&), const Matrix<T>& y) const;
        
    void unary(T (*fcn)(const T&), const Matrix<T>& x);
    void binary(T (*fcn)(const T&, const T&), const Matrix<T> &x, const Matrix<T> &y);
    void matrix_matrix(T (*fcn)(const T&, const T&), const Matrix<T>& x, const Matrix<T>& y);
    void matrix_scalar(T (*fcn)(const T&, const T&), const Matrix<T>& x, const T& y);
    void scalar_matrix(T (*fcn)(const T&, const T&), const T& x, const Matrix<T>& y);
#endif

    //@{
    /// Elementary operations -- Python naming
    Matrix<T> __add__(const Matrix<T> &y) const;
    Matrix<T> __sub__(const Matrix<T> &y) const;
    Matrix<T> __mul__(const Matrix<T> &y) const;
    Matrix<T> __div__(const Matrix<T> &y) const;
    Matrix<T> __pow__(const Matrix<T>& y) const;
    //@}

#ifndef SWIG
    /// Addition
    friend Matrix<T> operator+(const Matrix<T> &x, const Matrix<T> &y){ return x.__add__(y); }
    
    /// Subtraction
    friend Matrix<T> operator-(const Matrix<T> &x, const Matrix<T> &y){ return x.__sub__(y); }
    
    /// Elementwise multiplication
    friend Matrix<T> operator*(const Matrix<T> &x, const Matrix<T> &y){ return x.__mul__(y); }

    /// Elementwise division
    friend Matrix<T> operator/(const Matrix<T> &x, const Matrix<T> &y){ return x.__div__(y); }

    /// In-place addition
    Matrix<T>& operator+=(const Matrix<T> &y){return *this = this->__add__(y);}

    /// In-place subtraction
    Matrix<T>& operator-=(const Matrix<T> &y){return *this = this->__sub__(y);}

    /// In-place elementwise multiplication
    Matrix<T>& operator*=(const Matrix<T> &y){return *this = this->__mul__(y);}

    /// In-place elementwise division
    Matrix<T>& operator/=(const Matrix<T> &y){return *this = this->__div__(y);}
#endif // SWIG
    //@{
    /// Python operator overloading
    Matrix<T> __pow__ (const T& b) const{ return __pow__(Matrix<T>(b));}
    Matrix<T> __rpow__(const T& b) const{ return Matrix<T>(b).__pow__(*this);}
    Matrix<T> __add__ (const T& b) const{ return *this + b;}
    Matrix<T> __radd__(const T& b) const{ return b + *this;}
    Matrix<T> __sub__ (const T& b) const{ return *this - b;}
    Matrix<T> __rsub__(const T& b) const{ return b - *this;}
    Matrix<T> __mul__ (const T& b) const{ return *this * b;}
    Matrix<T> __rmul__(const T& b) const{ return b * *this;}
    Matrix<T> __div__ (const T& b) const{ return *this / b;}
    Matrix<T> __rdiv__(const T& b) const{ return b / *this;}
    //@}
    
    //@{
    /// Operations defined in the standard namespace for unambigous access and Numpy compatibility
    Matrix<T> sin() const;
    Matrix<T> cos() const;
    Matrix<T> tan() const;
    Matrix<T> arcsin() const;
    Matrix<T> arccos() const;
    Matrix<T> arctan() const;
    Matrix<T> exp() const;
    Matrix<T> log() const;
    Matrix<T> sqrt() const;
    Matrix<T> floor() const;
    Matrix<T> ceil() const;
    Matrix<T> fabs() const;
    Matrix<T> fmin(const Matrix<T>& y) const;
    Matrix<T> fmax(const Matrix<T>& y) const;
    //@}
    
    //@{
    /// Printing
#ifndef SWIG
    virtual void print(std::ostream &stream=std::cout) const; // print print description
    virtual void repr(std::ostream &stream=std::cout) const; // print representation
#endif
    std::string __repr__() { return getRepresentation(); } // python default print
    void printScalar(std::ostream &stream=std::cout) const; // print scalar
    void printVector(std::ostream &stream=std::cout) const; // print one row vector-style
    void printMatrix(std::ostream &stream=std::cout) const; // print one row, matrix-style
    void printSparse(std::ostream &stream=std::cout) const; // print sparse matrix style
    void printDense(std::ostream &stream=std::cout) const; // Print dense matrix stype
    //@}

    // Get the sparsity pattern
    const std::vector<int>& col() const;
    const std::vector<int>& rowind() const;
    int col(int el) const;
    int rowind(int row) const;
    void clear();
    void resize(int n, int m);
    void reserve(int nnz);
    void reserve(int nnz, int nrow);
    
    /** \brief Erase a submatrix
    Erase rows and/or columns of a matrix */
    void erase(const std::vector<int>& ii, const std::vector<int>& jj);
    
    /** \brief Enlarge matrix
    Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros */
    void enlarge(int nrow, int ncol, const std::vector<int>& ii, const std::vector<int>& jj);
    
    /// Access the non-zero elements
    std::vector<T>& data();
    
    /// Const access the non-zero elements
    const std::vector<T>& data() const;
    
    /// Const access the sparsity
    const CRSSparsity& sparsity() const;
    
    /// Access the sparsity, make a copy if there are multiple references to it
    CRSSparsity& sparsityRef();
    
    /** \brief  Set the non-zero elements, scalar */
    void set(T val, Sparsity sp=SPARSE);
    
    /** \brief  Get the non-zero elements, scalar */
    void get(T& val, Sparsity sp=SPARSE) const;

    /** \brief  Set the non-zero elements, vector */
    void set(const std::vector<T>& val, Sparsity sp=SPARSE);

    /** \brief  Get the non-zero elements, vector */
    void get(std::vector<T>& val, Sparsity sp=SPARSE) const;

    /** \brief  Set the non-zero elements, Matrix */
    void set(const Matrix<T>& val, Sparsity sp=SPARSE);

    /** \brief  Get the non-zero elements, Matrix */
    void get(Matrix<T>& val, Sparsity sp=SPARSE) const;

#ifdef SWIG
    %rename(get) getArray;
    %rename(set) setArray;
#endif

    /** \brief  Get the non-zero elements, array */
    void getArray(T* val, int len, Sparsity sp=SPARSE) const;

    /** \brief  Set the non-zero elements, array */
    void setArray(const T* val, int len, Sparsity sp=SPARSE);

    /** \brief  Get the non-zero elements, array, sparse and correct length */
    void getArray(T* val) const;

    /** \brief  Set the non-zero elements, array, sparse and correct length */
    void setArray(const T* val);
    
#ifndef SWIG
    /** \brief  Legacy - use getArray instead */
    void get(T* val, Sparsity sp=SPARSE) const;

    /** \brief  Legacy - use setArray instead */
    void set(const T* val, Sparsity sp=SPARSE);
#endif

    /** \brief  Get the result times a vector */
    void getTimesVector(const T *v, T *res) const;

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation 
    kl:    The number of subdiagonals in res 
    ku:    The number of superdiagonals in res 
    ldres: The leading dimension in res 
    res:   The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, T *res) const;

    // get the number if non-zeros for a given sparsity pattern
    int size(Sparsity sp) const;
    
  private:
    /// Sparsity of the matrix in a compressed row storage (CRS) format
    CRSSparsity sparsity_;
    
    /// Nonzero elements
    // std::vector<T> data_;
};

} // namespace CasADi

// Typedefs/template initializations
#ifndef SWIG
namespace CasADi{
  typedef Matrix<double> DMatrix;
  typedef std::vector<Matrix<double> > DMatrixVector;
  typedef std::vector< std::vector<Matrix<double> > > DMatrixVectorVector;
} // namespace CasADi
#else // SWIG
%template(DMatrix)             CasADi::Matrix<double>;
%template(DMatrixVector)       std::vector<CasADi::Matrix<double> > ;
%template(DMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<double> > > ;
#endif // SWIG

// The following functions must be placed in the standard namespace so that the old ones are not shadowed when CasADi namespace is used
#ifndef SWIG
namespace std{

  template<class T>
  CasADi::Matrix<T> sin(const CasADi::Matrix<T>& x){ return x.sin(); }

  template<class T>
  CasADi::Matrix<T> cos(const CasADi::Matrix<T>& x){ return x.cos(); }

  template<class T>
  CasADi::Matrix<T> tan(const CasADi::Matrix<T>& x){ return x.tan(); }

  template<class T>
  CasADi::Matrix<T> asin(const CasADi::Matrix<T>& x){ return x.arcsin(); }

  template<class T>
  CasADi::Matrix<T> acos(const CasADi::Matrix<T>& x){ return x.arccos(); }

  template<class T>
  CasADi::Matrix<T> atan(const CasADi::Matrix<T>& x){ return x.arctan(); }

  template<class T>
  CasADi::Matrix<T> exp(const CasADi::Matrix<T>& x){ return x.exp();}

  template<class T>
  CasADi::Matrix<T> log(const CasADi::Matrix<T>& x){ return x.log(); }

  template<class T>
  CasADi::Matrix<T> sqrt(const CasADi::Matrix<T>& x){ return x.sqrt();}

  template<class T>
  CasADi::Matrix<T> floor(const CasADi::Matrix<T>& x){ return x.floor();}

  template<class T>
  CasADi::Matrix<T> ceil(const CasADi::Matrix<T>& x){ return x.ceil();}

  template<class T>
  CasADi::Matrix<T> fabs(const CasADi::Matrix<T>& x){return x.fabs();}

  template<class T>
  CasADi::Matrix<T> fmin(const CasADi::Matrix<T>& x, const CasADi::Matrix<T>& y){ return x.fmin(y);}

  template<class T>
  CasADi::Matrix<T> fmax(const CasADi::Matrix<T>& x, const CasADi::Matrix<T>& y){ return x.fmax(y);}

  template<class T>
  CasADi::Matrix<T> pow(const CasADi::Matrix<T>& x, const CasADi::Matrix<T>& y){ return x.__pow__(y);}
  
} // namespace std
#endif // SWIG

#ifndef SWIG
namespace CasADi{
// Implementations

template<class T>
const T Matrix<T>::getElement(int i, int j) const{
  int ind = sparsity().getNZ(i,j);
  if(ind==-1)
    return 0;
  else
    return at(ind);
}

template<class T>
void Matrix<T>::setElement(int i, int j, const T& el){
  getElementRef(i,j) = el;
}

template<class T>
T& Matrix<T>::getElementRef(int i, int j){
  int oldsize = sparsity().size();
  int ind = sparsityRef().getNZ(i,j);
  if(oldsize != sparsity().size())
    data().insert(begin()+ind,0);
  return at(ind);
}

template<class T>
Matrix<T> Matrix<T>::getSub(const std::vector<int>& ii, const std::vector<int>& jj) const{
  // Nonzero mapping from submatrix to full
  std::vector<int> mapping;
  
  // Get the sparsity pattern
  CRSSparsity sp = sparsity().getSub(ii,jj,mapping);

  // Create return object
  Matrix<T> ret(sp);
  
  // Copy nonzeros
  for(int k=0; k<mapping.size(); ++k)
    ret[k] = (*this)[mapping[k]];
  
  // Return (RVO)
  return ret;
}

template<class T>
void Matrix<T>::setSub(const std::vector<int>& ii, const std::vector<int>& jj, const Matrix<T>& el){
  casadi_assert_message(el.numel()==1 || (ii.size() == el.size1() && jj.size() == el.size2()),"Dimension mismatch.");
  
  // If m is scalar
  if(el.numel() != ii.size() * jj.size()){
    setSub(ii,jj,Matrix<T>(ii.size(),jj.size(),el(0,0)));
    return;
  }

  if(dense() && el.dense()){
    // Dense mode
    for(int i=0; i<ii.size(); ++i) {
      for(int j=0; j<jj.size(); ++j) {
        (*this)[ii[i]*size2() + jj[j]]=el[i*el.size2()+j];
      }
    }
  } else {
    // Sparse mode

    // Remove submatrix to be replaced
    erase(ii,jj);

    // Extend el to the same dimension as this
    Matrix<T> el_ext = el;
    el_ext.enlarge(size1(),size2(),ii,jj);

    // Unite the sparsity patterns
    *this = unite(*this,el_ext);
  }
}

template<class T>
Matrix<T> Matrix<T>::getNZ(const std::vector<int>& kk) const{
  Matrix<T> ret(kk.size(),1,0);
  for(int k=0; k<kk.size(); ++k)
    ret[k] = data()[kk[k]];
  
  return ret;
}

template<class T>
void Matrix<T>::setNZ(const std::vector<int>& kk, const Matrix<T>& m){
  for(int k=0; k<kk.size(); ++k)
    data()[kk[k]] = m[k];
}

template<class T>
int Matrix<T>::size() const{
  return std::vector<T>::size();
}

template<class T>
int Matrix<T>::sizeU() const{
  return sparsity_.sizeU();
}

template<class T>
int Matrix<T>::sizeL() const{
  return sparsity_.sizeL();
}

template<class T>
int Matrix<T>::size1() const{
  return sparsity_.size1();
}

template<class T>
int Matrix<T>::size2() const{
  return sparsity_.size2();
}

template<class T>
int Matrix<T>::numel() const{
  return size1()*size2();
}

        
template<class T>
void Matrix<T>::makeDense(int n, int m, const T& val){
  if(n*m != numel()){
    sparsity_ = CRSSparsity(n,m,true);
    data().clear();
    data().resize(n*m, val);
  }
}

template<class T>
bool Matrix<T>::empty() const{
  return numel()==0;
}

template<class T>
bool Matrix<T>::scalar() const{
  return size1()==1 && size2()==1;
}

template<class T>
bool Matrix<T>::vector() const{
  return size2()==1;
}

template<class T>
bool Matrix<T>::dense() const{
  return size()==numel();
}

template<class T>
Matrix<T>::Matrix(){
  sparsity_ = CRSSparsity(0,0,false);
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& m) : std::vector<T>(m.data()){
  sparsity_ = m.sparsity_;
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x) : std::vector<T>(x){
  sparsity_ = CRSSparsity(x.size(),1,true);
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x, int n, int m) : std::vector<T>(x){
  casadi_assert_message(x.size() == n*m, "dimension mismatch");
  sparsity_ = CRSSparsity(n,m,true);
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
  sparsity_ = m.sparsity_;
  static_cast<std::vector<T>&>(*this) = m;
  return *this;
}

template<class T>
Matrix<T>::Matrix(int n, int m){
  sparsity_ = CRSSparsity(n,m,false);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T& val){
  sparsity_ = CRSSparsity(n,m,true);
  data().resize(n*m, val);
}

template<class T>
void Matrix<T>::makeEmpty(int n, int m){
  sparsity_ = CRSSparsity(n,m,false);
  data().clear();
}

template<class T>
void Matrix<T>::printScalar(std::ostream &stream) const {
  casadi_assert_message(numel()==1, "not a scalar");
  stream << (*this)(0);
}
  
template<class T>
void Matrix<T>::printVector(std::ostream &stream) const {
  // print dimension and first element
  stream << "[";
  stream << (*this)(0);
  for(int i=1; i<size1(); i++){
    stream << "," << (*this)(i);
  }
  stream << "]";  
}

template<class T>
void Matrix<T>::printMatrix(std::ostream &stream) const{ 
 // print the first row
 stream << "[[";
  stream << (*this)(0,0);
  for(int j=1; j<size2(); j++){
    stream << "," << (*this)(0,j);
 }
 stream << "]";
 // print the rest of the rows
 for(int i=1; i<size1(); i++){  
    stream << ",["; 
    stream << (*this)(i,0);
    for(int j=1; j<size2(); j++){
      stream << "," << (*this)(i,j);
    }
    stream << "]";
  }
  stream << "]";  
}

template<class T>
void Matrix<T>::printDense(std::ostream &stream) const{
  stream << size1() << "-by-" << size2() << " matrix:" << std::endl;
  for(int i=0; i<size1(); ++i){
    if(i==0)
      stream << "[[";
    else
      stream << " [";
    int j=0;
    for(int el=rowind(i); el<rowind(i+1); ++el){
      // Print leading zeros
      for(;j<col(el); ++j)
        stream << "0,  ";
      
      // Print element
      stream << (*this)[el] << ",  ";
    }
    
    // Print trailing zeros
    for(;j<size2(); ++j)
      stream << "0,  ";
    
    // New row
    if(i==size1()-1)
      stream << "]]" << std::endl;
    else
      stream << "]" << std::endl;
  }
}

template<class T>
void Matrix<T>::printSparse(std::ostream &stream) const {
  stream << "Sparse " << size1() << "-by-" << size2() << " matrix with " << size() << " structural non-zeros:" << std::endl;
  for(int i=0; i<size1(); ++i)
    for(int el=rowind(i); el<rowind(i+1); ++el){
      int j=col(el);
      stream << "(" << i << "," << j << "): " << (*this)[el] << std::endl;
    }
}

template<class T>
void Matrix<T>::print(std::ostream &stream) const{
  if(dense()){
    printDense(stream);
  } else {
    printSparse(stream);
  }
}

template<class T>
void Matrix<T>::repr(std::ostream &stream) const{
  if (empty())
    stream << "[]";
  else if (numel()==1)
    printScalar(stream);
  else if (size2()==1)
    printVector(stream);
  else
    printMatrix(stream);
}

template<class T>
const std::vector<int>& Matrix<T>::col() const{
  return sparsity().col();
}

template<class T>
const std::vector<int>& Matrix<T>::rowind() const{
  return sparsity_.rowind();
}

template<class T>
int Matrix<T>::col(int el) const{
  return sparsity_.col(el);
}

template<class T>
int Matrix<T>::rowind(int row) const{
  return sparsity_.rowind(row);
}

template<class T>
void Matrix<T>::reserve(int nnz){
  reserve(nnz,size1());
}

template<class T>
void Matrix<T>::reserve(int nnz, int nrow){
  data().reserve(nnz);
  sparsity_.reserve(nnz,nrow);
}

template<class T>
void Matrix<T>::resize(int n_, int m_){
  sparsity_.resize(n_,m_);
}

template<class T>
void Matrix<T>::clear(){
  sparsity_ = CRSSparsity(0,0,false);
  data().clear();
}

template<class T>
Matrix<T>::Matrix(const Element<Matrix<T>,T>& val){
  sparsity_ = CRSSparsity(1,1,true);
  data().resize(1,T(val));
}

template<class T>
Matrix<T>::Matrix(double val){
  sparsity_ = CRSSparsity(1,1,true);
  data().resize(1,val);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind, const std::vector<T>& data) : std::vector<T>(data){
  sparsity_ = CRSSparsity(n,m,col,rowind);
  if(data.size() != sparsity_.size())
    std::vector<T>::resize(sparsity_.size());
}

template<class T>
Matrix<T>::Matrix(const CRSSparsity& sparsity, const T& val){
  sparsity_ = sparsity;
  data().resize(sparsity_.size(),val);
}


template<class T>
const T Matrix<T>::getitem(int i) const{
  return data().at(i);
}

template<class T>
const T Matrix<T>::getitem(const std::vector<int> &I) const{
  casadi_assert_message(I.size()==2,"Index vector must be two-dimensional");
  return getElement(I[0],I[1]);
}

template<class T>
const Matrix<T> Matrix<T>::getitem(const std::vector< std::vector<int> > &II) const{
  casadi_assert_message(II.size()==2,"Index vector must be two-dimensional");
  return (*this)(II[0],II[1]);
}

template<class T>
void Matrix<T>::setitem(int k, const T& el){ 
  data().at(k) = el;
}

template<class T>
void Matrix<T>::setZero(){
  setAll(0);
}

template<class T>
void Matrix<T>::setAll(const T& val){
  std::fill(begin(),end(),val);
}

template<class T>
void Matrix<T>::setitem(const std::vector<int> &I, const T&  el){ 
  casadi_assert_message(I.size()==2,"Index vector must be two-dimensional");
  getElementRef(I[0],I[1]) = el;
}

template<class T>
void Matrix<T>::setitem(const std::vector< std::vector<int> > &II, const Matrix<T>& m){
  casadi_assert_message(II.size()==2,"Index vector must be two-dimensional");
  setSub(II[0],II[1],m);
}

template<class T>
Matrix<T> Matrix<T>::unary(T (*fcn)(const T&)) const{
  Matrix<T> temp;
  temp.unary(fcn,*this);
  return temp;
}

template<class T>
void Matrix<T>::unary(T (*fcn)(const T&), const Matrix<T>& x){
  // First check the value of the zero-entries
  T fcn_0 = fcn(0);
  if(casadi_limits<T>::isZero(fcn_0)){
    // Copy the matrix, including the sparsity pattern
    *this = x;
    
    // Do the operation on all non-zero elements
    for(int el=0; el<size(); ++el)
      (*this)[el] = fcn(x[el]);
    
  } else {
    makeDense(x.size1(),x.size2(),fcn_0);
    for(int i=0; i<size1(); ++i){ // loop over rows
      for(int el=x.rowind(i); el<x.rowind(i+1); ++el){ // loop over non-zero elements
        int j = x.col(el);
        (*this)[j+i*size2()] = fcn(x[el]);
      }
    }
  }
}

template<class T>
Matrix<T> Matrix<T>::binary(T (*fcn)(const T&, const T&), const Matrix<T>& y) const{
  Matrix<T> temp;
  temp.binary(fcn,*this,y);
  return temp;
}

template<class T>
void Matrix<T>::binary(T (*fcn)(const T&, const T&), const Matrix<T> &x, const Matrix<T> &y){
  if(x.scalar())
    if(y.scalar())
      *this = fcn(x(0),y(0));
    else
      scalar_matrix(fcn,x(0),y);
  else if(y.scalar())
    matrix_scalar(fcn,x,y(0));
  else
    matrix_matrix(fcn,x,y);
}

template<class T>
void Matrix<T>::scalar_matrix(T (*fcn)(const T&, const T&), const T& x, const Matrix<T>& y){
  T fcn_x_0 = fcn(x,0);
  if(casadi_limits<T>::isZero(fcn_x_0))
    // Start with an empty matrix: all elements are added to the end!
    makeEmpty(y.size1(),y.size2());
  else
    // Start with an dense matrix: element access in constant time
    makeDense(y.size1(),y.size2(),fcn_x_0);
  
  for(int i=0; i<size1(); ++i){ // loop over rows
    for(int el=y.rowind(i); el<y.rowind(i+1); ++el){
      int j = y.col(el);
      getElementRef(i,j) = fcn(x,y[el]);
    }
  }
}

template<class T>
void Matrix<T>::matrix_scalar(T (*fcn)(const T&, const T&), const Matrix<T>& x, const T& y){
  T fcn_0_y = fcn(0,y);
  if(casadi_limits<T>::isZero(fcn_0_y))
    // Start with an empty matrix: all elements are added to the end!
    makeEmpty(x.size1(),x.size2());
  else
    // Start with an dense matrix: element access in constant time
    makeDense(x.size1(),x.size2(),fcn_0_y);
  
  for(int i=0; i<size1(); ++i){ // loop over rows
    for(int el=x.rowind(i); el<x.rowind(i+1); ++el){
      int j = x.col(el);
      getElementRef(i,j) = fcn(x[el],y);
    }
  }
}

template<class T>
void Matrix<T>::matrix_matrix(T (*fcn)(const T&, const T&), const Matrix<T>& x, const Matrix<T>& y){
if(x.size1() != y.size1() || x.size2() != y.size2()) throw CasadiException("matrix_matrix: dimension mismatch");
  // Make a deep copy if *this and x or y is the same object
  if(this == &x)
    *this = x;
  else if(this == &y)
    *this = y;

  T fcn_0_0 = fcn(0,0);
  if(casadi_limits<T>::isZero(fcn_0_0))
    // Start with an empty matrix: all elements are added to the end!
    makeEmpty(x.size1(),x.size2());
  else
    // Start with an dense matrix: element access in constant time
    makeDense(x.size1(),x.size2(),fcn_0_0);
 
  for(int i=0; i<size1(); ++i){ // loop over rows
    int el1 = x.rowind(i);
    int el2 = y.rowind(i);
    int k1 = x.rowind(i+1);
    int k2 = y.rowind(i+1);
    while(el1 < k1 || el2 < k2){
      int j1 = (el1 < k1) ? x.col(el1) : numel() ;
      int j2 = (el2 < k2) ? y.col(el2) : numel() ;
      
      if(j1==j2)
        getElementRef(i,j1) = fcn(x[el1++],y[el2++]); 
      else if(j1>j2)
        getElementRef(i,j2) = fcn(0,y[el2++]);
      else
        getElementRef(i,j1) = fcn(x[el1++],0);
      }
    }
}

template<class T>
Matrix<T> Matrix<T>::operator-() const{
  Matrix<T> temp;
  temp.unary(casadi_operators<T>::neg,*this);
  return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator+() const{
  return *this;
}

template<class T>
Matrix<T> Matrix<T>::__add__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::add,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::__sub__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::sub,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::__mul__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::mul,*this,y);
  return r;
}

template<class T>
Matrix<T> Matrix<T>::__div__(const Matrix<T> &y) const{
  Matrix<T> r;
  r.binary(casadi_operators<T>::div,*this,y);
  return r;
}

// template<class T>
// Matrix<T>& Matrix<T>::operator+=(const double &y){
//   Matrix<T> x = *this;
//   binary(casadi_operators<T>::add,x,y);
//   return *this;
// }
// 
// template<class T>
// Matrix<T>& Matrix<T>::operator-=(const double &y){
//   Matrix<T> x = *this;
//   binary(casadi_operators<T>::sub,x,y);
//   return *this;
// }
// 
// template<class T>
// Matrix<T>& Matrix<T>::operator*=(const double &y){
//   Matrix<T> x = *this;
//   binary(casadi_operators<T>::mul,x,y);
//   return *this;
// }
// 
// template<class T>
// Matrix<T>& Matrix<T>::operator/=(const double &y){
//   Matrix<T> x = *this;
//   binary(casadi_operators<T>::div,x,y);
//   return *this;
// }

template<class T>
const CRSSparsity& Matrix<T>::sparsity() const{
  return sparsity_;
}

template<class T>
CRSSparsity& Matrix<T>::sparsityRef(){
  sparsity_.makeUnique();
  return sparsity_;
}

template<class T>
void Matrix<T>::getTimesVector(const T *v, T *res) const{
  // copy the result
  for(int i=0; i<size1(); ++i){ // loop over rows
    res[i] = 0;
    for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
      int j=col(el);  // column

      // Multiply with the vector
      res[i] += v[j]*(*this)[el];
    }
  }
}

template<class T>
void Matrix<T>::getBand(int kl, int ku, int ldres, T *res) const{
  // delete the content of the matrix
  for(int j=0; j<size2(); ++j) // loop over columns
    for(int s=0; s<kl+ku+1; ++s) // loop over the subdiagonals
      res[s + ldres*j] = 0;
  
  // loop over rows
  for(int i=0; i<size1(); ++i){ 
    
    // loop over the non-zero elements
    for(int el=rowind(i); el<rowind(i+1); ++el){ 
      int j=col(el);  // column
      
      // Check if we have not yet inside the band
      if(j<i-kl) continue;

      // Check if we are already outside the band
      if(j>i+ku) break;

      // Get the subdiagonal
      int s = i - j + ku;

      // Store the element
      res[s + ldres*j] = (*this)[el];
    }
  }
}

template<class T>
void Matrix<T>::set(T val, Sparsity sp){
  setArray(&val,1,DENSE);
}
    
template<class T>
void Matrix<T>::get(T& val, Sparsity sp) const{
  getArray(&val,1,DENSE);
}

template<class T>
void Matrix<T>::set(const std::vector<T>& val, Sparsity sp){
  setArray(&val[0],val.size(),sp);
}

template<class T>
void Matrix<T>::get(std::vector<T>& val, Sparsity sp) const{
  getArray(&val[0],val.size(),sp);
}

template<class T>
void Matrix<T>::set(const Matrix<T>& val, Sparsity sp){
  set(val.data(),sp);
}

template<class T>
void Matrix<T>::get(Matrix<T>& val, Sparsity sp) const{
  get(val.data(),sp);
}

template<class T>
void Matrix<T>::set(const T* val, Sparsity sp){
  int len = sp==SPARSE ? size() : sp==DENSE ? numel() : sp==SPARSESYM ? sizeL() : -1;
  setArray(val,len,sp);
}

template<class T>
void Matrix<T>::get(T* val, Sparsity sp) const{
  int len = sp==SPARSE ? size() : sp==DENSE ? numel() : sp==SPARSESYM ? sizeL() : -1;
  getArray(val,len,sp);
}

template<class T>
void Matrix<T>::getArray(T* val, int len, Sparsity sp) const{
  const std::vector<T> &v = *this;
  if(sp==SPARSE || (sp==DENSE && numel()==v.size())){
    if (len!=size()) {
			std::stringstream s;
		  s << "Matrix<T>::setArray: number of non-zero elements is expected to be " << size() << ", but got " << len << " instead.";
      throw CasadiException(s.str());
    }
    copy(v.begin(),v.end(),val);
  } else if(sp==DENSE){
    if (len!=numel()) {
			std::stringstream s;
		  s << "Matrix<T>::setArray: number of elements is expected to be " << numel() << ", but got " << len << " instead.";
      throw CasadiException(s.str());
    }
    int k=0; // index of the result
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
        int j=col(el);  // column
        for(; k<i*size2()+j; ++k)
          val[k] = 0; // add zeros before the non-zero element

      // add the non-zero element
      val[k] = v[el];
      k++;
    }
    // add sparse zeros at the end of the matrix
    for(; k<numel(); ++k)
     val[k] = 0;
  } else if(sp==SPARSESYM){
    // copy to the result vector
    int nz = 0;
    for(int row=0; row<size1(); ++row){
      // Loop over the elements in the row
      for(int el=rowind(row); el<rowind(row+1); ++el){ // loop over the non-zero elements
        if(col(el) > row) break; // break inner loop (only lower triangular part is used)
        val[nz++] = v[el];
      }
    }
  } else {
    throw CasadiException("Matrix<T>::getArray: not SPARSE or DENSE");
  }
}

template<class T>
void Matrix<T>::setArray(const T* val, int len, Sparsity sp){
  std::vector<T> &v = *this;
  if(sp==SPARSE || (sp==DENSE && numel()==size())){
    if (len!=size()) {
			std::stringstream s;
		  s << "Matrix<T>::setArray: number of non-zero elements is expected to be " << size() << ", but got " << len << " instead.";
      throw CasadiException(s.str());
    }
    copy(val,val+len,v.begin());
  } else if(sp==DENSE){
    if (len!=numel()) {
			std::stringstream s;
		  s << "Matrix<T>::setArray: number of elements is expected to be " << numel() << ", but got " << len << " instead.";
      throw CasadiException(s.str());
    }
    for(int i=0; i<size1(); ++i) // loop over rows
      for(int el=rowind(i); el<rowind(i+1); ++el){ // loop over the non-zero elements
        // column
        int j=col(el);
        
        // Set the element
        v[el] = val[i*size2()+j];
    }
  } else {
    throw CasadiException("Matrix<T>::setArray: not SPARSE or DENSE");
  }
}

template<class T>
void Matrix<T>::getArray(T* val) const{
  getArray(val,size(),SPARSE);
}

template<class T>
void Matrix<T>::setArray(const T* val){
  setArray(val,size(),SPARSE);
}






template<class T>
int Matrix<T>::size(Sparsity sp) const{
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

template<class T>
Matrix<T> Matrix<T>::__pow__(const Matrix<T>& y) const{
  return binary(CasADi::casadi_operators<T>::pow,y);
}

template<class T>
Matrix<T> Matrix<T>::sin() const{
  return unary(CasADi::casadi_operators<T>::sin);
}

template<class T>
Matrix<T> Matrix<T>::cos() const{
  return unary(CasADi::casadi_operators<T>::cos);
}

template<class T>
Matrix<T> Matrix<T>::tan() const{
  return unary(CasADi::casadi_operators<T>::tan);
}

template<class T>
Matrix<T> Matrix<T>::arcsin() const{
  return unary(CasADi::casadi_operators<T>::asin);
}

template<class T>
Matrix<T> Matrix<T>::arccos() const{
  return unary(CasADi::casadi_operators<T>::acos);
}

template<class T>
Matrix<T> Matrix<T>::arctan() const{
  return unary(CasADi::casadi_operators<T>::atan);
}

template<class T>
Matrix<T> Matrix<T>::exp() const{
  return unary(CasADi::casadi_operators<T>::exp);
}

template<class T>
Matrix<T> Matrix<T>::log() const{
  return unary(CasADi::casadi_operators<T>::log);
}

template<class T>
Matrix<T> Matrix<T>::sqrt() const{
  return unary(CasADi::casadi_operators<T>::sqrt);
}

template<class T>
Matrix<T> Matrix<T>::floor() const{
  return unary(CasADi::casadi_operators<T>::floor);
}

template<class T>
Matrix<T> Matrix<T>::ceil() const{
  return unary(CasADi::casadi_operators<T>::ceil);
}

template<class T>
Matrix<T> Matrix<T>::fabs() const{
  return unary(CasADi::casadi_operators<T>::fabs);
}

template<class T>
Matrix<T> Matrix<T>::fmin(const Matrix<T>& y) const{
  return binary(CasADi::casadi_operators<T>::fmin, y);
}

template<class T>
Matrix<T> Matrix<T>::fmax(const Matrix<T>& y) const{
  return binary(CasADi::casadi_operators<T>::fmax, y);
}

template<class T>
std::vector<T>& Matrix<T>::data(){
  return *this;  
}
    
template<class T>
const std::vector<T>& Matrix<T>::data() const{
  return *this;  
}

template<class T>
void Matrix<T>::erase(const std::vector<int>& ii, const std::vector<int>& jj){
  // Erase from sparsity pattern
  std::vector<int> mapping = sparsityRef().erase(ii,jj);
  
  // Update non-zero entries
  for(int k=0; k<mapping.size(); ++k)
    (*this)[k] = (*this)[mapping[k]];
    
  // Truncate nonzero vector
  data().resize(mapping.size());
}

template<class T>
void Matrix<T>::enlarge(int nrow, int ncol, const std::vector<int>& ii, const std::vector<int>& jj){
  sparsityRef().enlarge(nrow,ncol,ii,jj);
}


} // namespace CasADi
#endif // SWIG



#endif // MATRIX_HPP

