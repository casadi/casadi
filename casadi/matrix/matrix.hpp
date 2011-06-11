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
#include "slice.hpp"
#include "element.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "crs_sparsity.hpp"

#if defined(_WIN32) || defined(_WIN64)
#pragma warning (disable:4996)
#pragma warning (disable:4018)
#endif

namespace CasADi{
  
  /** Sparsity format for getting and setting inputs and outputs */
  enum Sparsity{SPARSE,SPARSESYM,DENSE,DENSESYM};

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
class Matrix : public PrintableObject{

  public:
    /** \brief  constructors */
    /// empty 0-by-0 matrix constructor
    Matrix();
    
    /// Copy constructor
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
    Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind, const std::vector<T>& d=std::vector<T>());

    /// dense matrix constructor with data given as vector of vectors
    explicit Matrix(const std::vector< std::vector<T> >& m);
    
    /// sparse matrix with a given sparsity
    explicit Matrix(const CRSSparsity& sparsity, const T& val=0);
    
    /// sparse matrix with a given sparsity and non-zero elements.
    Matrix(const CRSSparsity& sparsity, const std::vector<T>& d);
    
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

    /** \brief  Create a matrix from a matrix with a different type of matrix entries (assuming that the scalar conversion is valid) */
    template<typename A>
    Matrix(const Matrix<A>& x){
      sparsity_ = x.sparsity();
      data().resize(x.size());
      copy(x.begin(),x.end(),begin());
    }

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
      data().resize(numel());
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

#endif // SWIG

    /// Get nonzeros
    Matrix<T> getNZ(const std::vector<int>& kk) const;
    
    /** \brief set nonzeros
    m.size() must equal kk.size() unless m is scalar
    */
    void setNZ(const std::vector<int>& kk, const Matrix<T>& m);

    void setNZ(const std::vector<int>& kk, const T& m);
    
#ifndef SWIG 

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
    const Matrix<T> operator()(const Slice& i, A j) const{ return operator()(i.getAll(size1()),j);}
    
    /// Get all rows
    template<class A>
    SubMatrix<Matrix<T> > operator()(const Slice& i, A j){ return operator()(i.getAll(size1()),j);}
  
    /// Access all columns
    template<class A>
    const Matrix<T> operator()(A i, const Slice& j) const{ return operator()(i,j.getAll(size2()));}

    /// Get all columns
    template<class A>
    SubMatrix<Matrix<T> > operator()(A i, const Slice& j){ return operator()(i,j.getAll(size2()));}

    /// Get all rows and columns
    const Matrix<T> operator()(const Slice& i, const Slice& j) const{ return operator()(i.getAll(size1()),j.getAll(size2()));}

    /// Access all rows and columns
    SubMatrix<Matrix<T> > operator()(const Slice& i, const Slice& j){ return operator()(i.getAll(size1()),j.getAll(size2()));}

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
    const T __getitem__(int i) const;
    
    /// Get several non-zero entries
    const Matrix<T> __getitem__(const Slice& k) const;
    
    /// Get several non-zero entries
    const Matrix<T> __getitem__(const std::vector<int>& k) const;

    /// Python: get a matrix entry
    const T __getitem__(int i, int j) const;
    
    /// Python: get a submatrix
    const Matrix<T> __getitem__(const std::vector<Slice> & ij) const;
    
    /// Python: get a submatrix
    const Matrix<T> __getitem__(const std::pair<std::vector<int>,std::vector<int> > & ij) const;
      
    /// Python: set a non-zero entry
    void __setitem__(int k, const T& el);
    
    /// Python: set several non-zero entries
    void __setitem__(const Slice& k, const Matrix<T>& m);
 
    /// Python: set several non-zero entries
    void __setitem__(const std::vector<int>& k, const Matrix<T>& m);
       
    /// Python: set a matrix entry
    void __setitem__(int i, int j, const T&  el);

    /// Python: set a submatrix
    void __setitem__(const std::vector<Slice> & ij, const Matrix<T>& m);

    /// Python: set a submatrix
    void __setitem__(const std::pair<std::vector<int>, std::vector<int> > &j, const Matrix<T>& m);
    
    #ifdef SWIGOCTAVE
    
    /// Octave: get a non-zero
    T __paren__(int k) const{ return at(k-1);}
    Matrix<T> __paren__(const IndexList &k) const{
      return (*this)[k.getAll(size())];
    }
    Matrix<T> __paren__(const Slice &k) const{ 
      return (*this)[k.getAll(size())];
    }
    
    /// Octave: get a matrix element
    T __paren__(int i, int j) const{ return (*this)(i-1,j-1);}
    Matrix<T> __paren__(const IndexList &i, const IndexList &j) const{ 
      return (*this)(i.getAll(size1()),j.getAll(size2()));
    }
    Matrix<T> __paren__(const Slice &i, const Slice &j) const{ 
      return (*this)(i.getAll(size1()),j.getAll(size2()));
    }
    
    /// Octave: set a non-zero
    void __paren_asgn__(int k, const Matrix<T>& m){ at(k-1) = m(0,0);}
    void __paren_asgn__(const IndexList &k, const Matrix<T>& m){
      (*this)[k.getAll(size())] = m;
    }
    void __paren_asgn__(const Slice &k, const Matrix<T>& m){
      (*this)[k.getAll(size())] = m;
    }
    
    /// Octave: set a matrix element
    void __paren_asgn__(int i, int j, const Matrix<T>& m){ (*this)(i-1,j-1) = m(0,0);}
    void __paren_asgn__(const IndexList &i, const IndexList &j, const Matrix<T>& m){
      (*this)(i.getAll(size1()),j.getAll(size2())) = m;
    }
    
    void __paren_asgn__(const Slice &i, const Slice &j, const Matrix<T>& m){
      (*this)(i.getAll(size1()),j.getAll(size2())) = m;
    }
    
    #endif // SWIGOCTAVE
    
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
    /// Elementwise operations -- Octave/Python naming
    Matrix<T> __add__(const Matrix<T> &y) const;
    Matrix<T> __sub__(const Matrix<T> &y) const;
    Matrix<T> __mul__(const Matrix<T> &y) const;
    Matrix<T> __div__(const Matrix<T> &y) const;
    Matrix<T> __pow__(const Matrix<T> &y) const;
    Matrix<T> __radd__(const Matrix<T> &y) const {return y.__add__(*this);}
    Matrix<T> __rsub__(const Matrix<T> &y) const {return y.__sub__(*this);}
    Matrix<T> __rmul__(const Matrix<T> &y) const {return y.__mul__(*this);}
    Matrix<T> __rdiv__(const Matrix<T> &y) const {return y.__div__(*this);}
    Matrix<T> __rpow__(const Matrix<T> &y) const {return y.__pow__(*this);}
    //@}
    
    /// Matrix product
    Matrix<T> prod(const Matrix<T> &y) const;

    /// Matrix transpose
    Matrix<T> trans() const;

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
    
    /** \brief Get string representation of dimensions.
    The representation is (nrow x ncol = numel | size)
    */
    std::string dimString() const;

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
    %rename(get) getStridedArray;
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
    
    /** \brief  Get the non-zero elements, strided array */
    void getStridedArray(T* val, int len, int stride1, int stride2, Sparsity sp=SPARSE) const;
    
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
    std::vector<T> data_;
};

} // namespace CasADi

// Typedefs/template initializations
namespace CasADi{
  typedef Matrix<double> DMatrix;
  typedef std::vector<Matrix<double> > DMatrixVector;
  typedef std::vector< std::vector<Matrix<double> > > DMatrixVectorVector;
} // namespace CasADi

#ifdef SWIG // SWIG
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
void Matrix<T>::setNZ(const std::vector<int>& kk, const T& m){
  for(int k=0; k<kk.size(); ++k)
    data()[kk[k]] = m;
}

template<class T>
void Matrix<T>::setNZ(const std::vector<int>& kk, const Matrix<T>& m){
  if (m.size()==1 && m.numel()==1) {
    // Allow scalar assignment:
    // m[:2]=3
    for(int k=0; k<kk.size(); ++k)
      data()[kk[k]] = m[0];
    return;
  }
  if (kk.size()!=m.size()) {
    std::stringstream ss;
    ss << "Matrix<T>::setNZ: length of non-zero indices (" << kk.size() << ") " << std::endl;
    ss << "must match size of rhs (" << m.size() << ")." << std::endl;
    throw CasadiException(ss.str());
  }
  for(int k=0; k<kk.size(); ++k)
    data()[kk[k]] = m[k];
}

template<class T>
int Matrix<T>::size() const{
  return data().size();
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
  // Quick return if already dense
  if(n*m == size())
    return;
  
  if(size1()!=n || size2()!=m){
    // Also resize
    sparsity_ = CRSSparsity(n,m,true);
    std::fill(data().begin(),data().end(),val);
    data().resize(n*m, val);
  } else {
    // Create a new data vector
    data().resize(n*m,val);
    
    // Loop over the rows in reverse order
    for(int i=n-1; i>=0; --i){
      // Loop over nonzero elements in reverse order
      for(int el=rowind(i+1)-1; el>=rowind(i); --el){
        // Column
        int j = col(el);
        
        // Swap the old position with the new position
        if(el!=j+i*m){
          data()[j+i*m] = data()[el];
          data()[el] = val;
        }
      }
    }
      
    // Save the new sparsity pattern
    sparsity_ = CRSSparsity(n,m,true);
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
Matrix<T>::Matrix(const Matrix<T>& m){
  data_ = m.data_;
  sparsity_ = m.sparsity_;
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x) : data_(x){
  sparsity_ = CRSSparsity(x.size(),1,true);
}

template<class T>
Matrix<T>::Matrix(const std::vector<T>& x, int n, int m) : data_(x){
  casadi_assert_message(x.size() == n*m, "dimension mismatch");
  sparsity_ = CRSSparsity(n,m,true);
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
  sparsity_ = m.sparsity_;
  data_ = m.data_;
  return *this;
}

template<class T>
Matrix<T>::Matrix(int n, int m){
  sparsity_ = CRSSparsity(n,m,false);
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T& val){
  sparsity_ = CRSSparsity(n,m,true);
  data_.resize(n*m, val);
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
      j++;
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
Matrix<T>::Matrix(int n, int m, const std::vector<int>& col, const std::vector<int>& rowind, const std::vector<T>& d) : data_(d){
  sparsity_ = CRSSparsity(n,m,col,rowind);
  if(data_.size() != sparsity_.size())
    data_.resize(sparsity_.size());
}

template<class T>
Matrix<T>::Matrix(const std::vector< std::vector<T> >& d){
  int n=d.size();
  int m=-1;
  for (int i=0;i<n;i++) {
    if (m==-1) {
      m = d[i].size();
    } else {
      if (m!=d[i].size()) {
        std::stringstream s;
		    s << "Matrix<T>::Matrix(const std::vector< std::vector<T> >& d): shape mismatch" << std::endl;
		    s << "Attempting to construct a matrix from a nested list." << std::endl;
		    s << "I got convinced that the desired size is ("<< n << " x " << m << " ), but now I encounter a vector of size (" << d[i].size() <<  " )" << std::endl;
		    throw CasadiException(s.str());
      }
    }
  }
  sparsity_ = CRSSparsity(n,m,true);
  
  data().resize(n*m);

  for (int i=0;i<n;i++) {
    copy(d[i].begin(),d[i].end(),begin()+i*m);
  }
  
}

template<class T>
Matrix<T>::Matrix(const CRSSparsity& sparsity, const T& val){
  sparsity_ = sparsity;
  data().resize(sparsity_.size(),val);
}

template<class T>
Matrix<T>::Matrix(const CRSSparsity& sparsity, const std::vector<T>& d) : data_(d) {
  casadi_assert_message(sparsity.size()==d.size(),"size mismatch");
  sparsity_ = sparsity;
}

template<class T>
const T Matrix<T>::__getitem__(int i) const{
  return data().at(i);
}

template<class T>
const T Matrix<T>::__getitem__(int i, int j) const{
  return getElement(i,j);
}

template<class T>
const Matrix<T> Matrix<T>::__getitem__(const std::vector<Slice> & ij) const{
  casadi_assert_message(ij.size()==2,"Can only do up to 2D slices");
  // Note: we do not use pair, because that is very difficult to typemap (because of the comma in the type description)
  return (*this)(ij[0].getAll(size1()),ij[1].getAll(size2()));
}

template<class T>
const Matrix<T> Matrix<T>::__getitem__(const std::pair<std::vector<int>,std::vector<int> > & ij) const{
  return (*this)(ij.first,ij.second);
}

template<class T>
const Matrix<T> Matrix<T>::__getitem__(const Slice& kk) const{
  return (*this)[kk.getAll(size())];
}

template<class T>
const Matrix<T> Matrix<T>::__getitem__(const std::vector<int>& kk) const{
  return (*this)[kk];
}

template<class T>
void Matrix<T>::__setitem__(int k, const T& el){ 
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
void Matrix<T>::__setitem__(const Slice& k, const Matrix<T>& m){
  (*this)[k.getAll(size())] = m;
}

template<class T>
void Matrix<T>::__setitem__(const std::vector<int>& k, const Matrix<T>& m){
  (*this)[k] = m;
}

template<class T>
void Matrix<T>::__setitem__(int i, int j, const T&  el){ 
  getElementRef(i,j) = el;
}

template<class T>
void Matrix<T>::__setitem__(const std::vector<Slice> & ij, const Matrix<T>& m){
  casadi_assert_message(ij.size()==2,"Can only do up to 2D slices");
  // Note: we do not use pair, because that is very difficult to typemap (because of the comma in the type description)
  setSub(ij[0].getAll(size1()),ij[1].getAll(size2()),m);
}

template<class T>
void Matrix<T>::__setitem__(const std::pair<std::vector<int>, std::vector<int> > & ij, const Matrix<T>& m){
  setSub(ij.first,ij.second,m);
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
  const std::vector<T> &v = data();
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

/**
Set stride to zero for unstrided acces
*/
template<class T>
void Matrix<T>::getStridedArray(T* val, int len, int stride1, int stride2, Sparsity sp) const{
  if (stride1==0 || stride2==0 || (stride2==1 && stride1==size2())) 
    return getArray(val, len, sp);
    
  const std::vector<T> &v = data();
  if(sp==SPARSE){
    throw CasadiException("Matrix<T>::getArray: strided SPARSE not implemented");
  } else if(sp==DENSE && numel()==v.size()) {
    for (int i=0;i<size1();i++) {
      for (int j=0;j<size2();j++) {
        val[i*stride1+j*stride2] = v[i*size2()+j];
      }
    }
  } else if(sp==DENSE){
    throw CasadiException("Matrix<T>::getArray: strided sparse DMatrix->dense not implemented");
  } else if(sp==SPARSESYM){
    throw CasadiException("Matrix<T>::getArray: strided SPARSESYM not implemented");
  } else {
    throw CasadiException("Matrix<T>::getArray: not SPARSE or DENSE");
  }

}

template<class T>
void Matrix<T>::setArray(const T* val, int len, Sparsity sp){
  std::vector<T> &v = data();
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
  return data_;  
}
    
template<class T>
const std::vector<T>& Matrix<T>::data() const{
  return data_;  
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

template<class T>
std::string Matrix<T>::dimString() const {
  std::stringstream ss;
  ss << "(" << size1() << " x " << size2() << " = " << numel() << " | " << size() << ")";
  return ss.str();
}

template<class T>
Matrix<T> Matrix<T>::prod(const Matrix<T> &y) const{
  // First factor
  const Matrix<T>& x = *this;
  
  if (x.size2() != y.size1()) {
    std::stringstream ss;
    ss << "Matrix<T>::prod: dimension mismatch. Attemping product of (" << x.size1() << " x " << x.size2() << ") " << std::endl;
    ss << "with (" << y.size1() << " x " << y.size2() << ") matrix." << std::endl;
    throw CasadiException(ss.str());
  }

  
  // Find the mapping corresponding to the transpose of y (no need to form the transpose explicitly)
  std::vector<int> y_trans_map;
  CRSSparsity y_trans_sparsity = y.sparsity().transpose(y_trans_map);

  // Create the sparsity pattern for the matrix-matrix product
  CRSSparsity spres = x.sparsity().patternProduct(y_trans_sparsity);
  
  // Create the return object
  Matrix<T> ret(spres);
  
  // Direct access to the arrays
  const std::vector<int> &r_col = ret.col();
  const std::vector<int> &r_rowind = ret.rowind();
  const std::vector<int> &x_col = x.col();
  const std::vector<int> &y_row = y_trans_sparsity.col();
  const std::vector<int> &x_rowind = x.rowind();
  const std::vector<int> &y_colind = y_trans_sparsity.rowind();

  // loop over the row of the resulting matrix)
  for(int i=0; i<ret.size1(); ++i){
    for(int el=r_rowind[i]; el<r_rowind[i+1]; ++el){ // loop over the non-zeros of the resulting matrix
      int j = r_col[el];
      int el1 = x_rowind[i];
      int el2 = y_colind[j];
      ret[el]=0;
      while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
        int j1 = x_col[el1];
        int i2 = y_row[el2];      
        if(j1==i2){
          ret[el] += x[el1++] * y[y_trans_map[el2++]];
        } else if(j1<i2) {
          el1++;
        } else {
          el2++;
        }
      }
    }
  }
  return ret;
}

template<class T>
Matrix<T> Matrix<T>::trans() const{
  // quick return if empty or scalar
  if(empty() || scalar()) return *this;

  // Create the new sparsity pattern and the mapping
  std::vector<int> mapping;
  CRSSparsity s = sparsity().transpose(mapping);

  // create the return matrix
  Matrix<T> ret(s);
  
  // Copy the content
  for(int i=0; i<mapping.size(); ++i)
    ret[i] = (*this)[mapping[i]];
  
  return ret;
}

} // namespace CasADi
#endif // SWIG



#endif // MATRIX_HPP

