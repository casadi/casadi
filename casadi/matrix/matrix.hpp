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
#include <typeinfo>
#include "../casadi_exception.hpp"
#include "../printable_object.hpp"
#include "../casadi_limits.hpp"
#include "../casadi_operators.hpp"
#include "../stl_vector_tools.hpp"
#include "slice.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "crs_sparsity.hpp"
#include "../casadi_calculus.hpp"

namespace CasADi{
  
  /** Sparsity format for getting and setting inputs and outputs */
  enum Sparsity{SPARSE,SPARSESYM,DENSE,DENSESYM};
  
  //@{
  /** \brief Get typename */
  template <typename T> inline const char* typeName() { return typeid(T).name(); }
  template<> inline const char* typeName<double>() { return "double"; }
  template<> inline const char* typeName<float>() { return "float"; }
  template<> inline const char* typeName<int>() { return "int"; }
  template<> inline const char* typeName<long>() { return "long"; }
  //@}

  /** \brief General sparse matrix class
  General sparse matrix class that is designed with the idea that "everything is a matrix", that is, also scalars and vectors.\n
  This philosophy makes it easy to use and to interface in particularily with Python and Matlab/Octave.\n
  
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
    
    /** \brief Check if the dimensions and rowind,col vectors are compatible.
    * \param complete  set to true to also check elementwise
    * throws an error as possible result
    */
    void sanityCheck(bool complete=false) const;
    
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

    /// Convert to scalar type
    const T toScalar() const;
//    operator const T() const;
    
    /// Scalar type
    typedef T ScalarType;
    
#ifndef SWIG
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

    #ifndef SWIG
    /// Get the shape
    std::pair<int,int> shape() const;
    #endif

    //@{
    /// Check type of matrix
    bool empty() const; // is the matrix empty
    bool scalar() const; // is the matrix scalar
    bool vector() const; // is the matrix a vector
    bool dense() const; // is the matrix dense
    //@}


    #ifndef SWIG
    /// Get a non-zero element
    const T& at(int k) const{ if (k<0) k+=size(); return data().at(k); }
    
    /// Access a non-zero element
    T& at(int k){if (k<0) k+=size(); return data().at(k); }
    #endif // SWIG
    
    #ifdef SWIG
    /// Access a non-zero element
    T at(int k){if (k<0) k+=size(); return data().at(k); }
    #endif // SWIG
    
    #ifndef SWIG
    /// get an element
    const T& elem(int i, int j=0) const;
    
    /// get a reference to an element
    T& elem(int i, int j=0);
    #endif // SWIG
    
    #ifdef SWIG
    /// Access a non-zero element
    T elem(int i, int j=0) { return elem(i,j);}
    #endif // SWIG

    /// get an element, do not allocate
    const T getElement(int i, int j=0) const{ return elem(i,j);}
    
    /// Returns true if the matrix has a non-zero at location i,j
    bool hasNZ(int i, int j) const { return sparsity().hasNZ(i,j); }

    //@{
    /// Get a submatrix
    const Matrix<T> getSub(int i, int j) const;
    const Matrix<T> getSub(int i, const std::vector<int>& j) const{ return getSub(std::vector<int>(1,i),j);}
    const Matrix<T> getSub(const std::vector<int>& i, int j) const{ return getSub(i,std::vector<int>(1,j));}
    const Matrix<T> getSub(const std::vector<int>& i, const std::vector<int>& j) const;
    const Matrix<T> getSub(const Slice& i, const Slice& j) const{ return getSub(i.getAll(size1()),j.getAll(size2()));}
    const Matrix<T> getSub(int i, const Slice& j) const{ return getSub(std::vector<int>(1,i),j.getAll(size2()));}
    const Matrix<T> getSub(const Slice& i, int j) const{ return getSub(i.getAll(size1()),std::vector<int>(1,j));}
    const Matrix<T> getSub(const std::vector<int>& i, const Matrix<int>& k) const;
    const Matrix<T> getSub(const Matrix<int>& k, const std::vector<int>& j) const;
    const Matrix<T> getSub(const Slice& i, const Matrix<int>& k) const {return getSub(i.getAll(size1()),k);}
    const Matrix<T> getSub(const Matrix<int>& k, const Slice& j) const {return getSub(k,j.getAll(size2()));}
    const Matrix<T> getSub(const Matrix<int>& i, const Matrix<int>& j) const;
    const Matrix<T> getSub(const CRSSparsity& sp, int dummy = 0) const;
    //@}

    //@{
    /// Set a submatrix
    void setSub(int i, int j, const Matrix<T>& m);
    void setSub(int i, const std::vector<int>& j, const Matrix<T>& m){ setSub(std::vector<int>(1,i),j,m);}
    void setSub(const std::vector<int>& i, int j, const Matrix<T>& m){ setSub(i,std::vector<int>(1,j),m);}
    void setSub(const std::vector<int>& i, const std::vector<int>& j, const Matrix<T>& m);
    void setSub(const Slice& i, const Slice& j, const Matrix<T>& m){ setSub(i.getAll(size1()),j.getAll(size2()),m);}
    void setSub(const std::vector<int>& i, const Matrix<int>& k, const Matrix<T>& m);
    void setSub(const Matrix<int>& k, const std::vector<int>& j, const Matrix<T>& m);
    void setSub(const Slice& i, const Matrix<int>& k, const Matrix<T>& m) {return setSub(i.getAll(size1()),k,m);}
    void setSub(const Matrix<int>& k, const Slice& j, const Matrix<T>& m) {return setSub(k,j.getAll(size2()),m);}
    void setSub(const Matrix<int>& i, const Matrix<int>& j, const Matrix<T>& m);
    void setSub(const CRSSparsity& sp, int dummy, const Matrix<T>& m);
    //@}

    //@{
    /// Get a set of nonzeros
    const Matrix<T> getNZ(int k) const;
    const Matrix<T> getNZ(const std::vector<int>& k) const;
    const Matrix<T> getNZ(const Slice& k) const{ return getNZ(k.getAll(size()));}
    const Matrix<T> getNZ(const Matrix<int>& k) const;
    //@}
    
    //@{
    /// Set a set of nonzeros
    void setNZ(int k, const Matrix<T>& m);
    void setNZ(const std::vector<int>& k, const Matrix<T>& m);
    void setNZ(const Slice& k, const Matrix<T>& m){ setNZ(k.getAll(size()),m);}
    void setNZ(const Matrix<int>& k, const Matrix<T>& m);
    //@}

    /// Append a matrix to the end
    void append(const Matrix<T>& y);

#ifndef SWIG 
    /** \brief  Get vector nonzero or slice of nonzeros */
    template<typename K>
    const Matrix<T> operator[](const K& k) const{ return getNZ(k); }

    /** \brief  Access vector nonzero or slice of nonzeros */
    template<typename K>
    NonZeros<Matrix<T>,K> operator[](const K& k){ return NonZeros<Matrix<T>,K>(*this,k); }

    /** \brief  Get vector element or slice */
    template<typename I>
    const Matrix<T> operator()(const I& i) const{ return getSub(i,0);}

    /** \brief  Get Sparsity slice */
    const Matrix<T> operator()(const CRSSparsity& sp) const{ return getSub(sp); }
    
    /** \brief  Get Matrix element or slice */
    template<typename I, typename J>
    const Matrix<T> operator()(const I& i, const J& j) const{ return getSub(i,j); }

    /** \brief  Access vector element or slice */
    template<typename I>
    SubMatrix<Matrix<T>,I,int> operator()(const I& i){ return SubMatrix<Matrix<T>,I,int>(*this,i,0); }
       
    /** \brief  Access Matrix element or slice */
    template<typename I, typename J>
    SubMatrix<Matrix<T>,I,J> operator()(const I& i, const J& j){ return SubMatrix<Matrix<T>,I,J>(*this,i,j); }
    
#endif // SWIG
    
    //@{
    /// Indexing for interfaced languages
    /// get a non-zero
    const Matrix<T> indexed_one_based(int k) const{ return operator[](k-1);}
    const Matrix<T> indexed_zero_based(int k) const{ return operator[](k);}
    const Matrix<T> indexed_one_based(const Matrix<int>& k) const{ return operator[](k-1);}
    const Matrix<T> indexed_zero_based(const Matrix<int>& k) const{ return operator[](k);}
    const Matrix<T> indexed(const Slice &k) const{ return (*this)[k];}
    const Matrix<T> indexed(const IndexList &k) const{
      return (*this)[k.getAll(size())];
    }
    
    /// get a matrix element
    const Matrix<T> indexed_one_based(int i, int j) const{ return operator()(i-1,j-1);}
    const Matrix<T> indexed_zero_based(int i, int j) const{ return operator()(i,j);}
    const Matrix<T> indexed(const Slice &i, const Slice &j) const{ return (*this)(i,j); }
    const Matrix<T> indexed(const IndexList &i, const IndexList &j) const{ 
      return (*this)(i.getAll(size1()),j.getAll(size2()));
    }
    const Matrix<T> indexed(const Slice &i, const Matrix<int>& k) const{ return (*this)(i,k); }
    const Matrix<T> indexed(const IndexList &i, const Matrix<int>& k) const{ 
      return (*this)(i.getAll(size1()),k);
    }
    const Matrix<T> indexed(const Matrix<int>& k, const Slice &j) const{ return (*this)(k,j); }
    const Matrix<T> indexed(const Matrix<int>& k, const IndexList &j) const{ 
      return (*this)(k,j.getAll(size2()));
    }
    const Matrix<T> indexed(const Matrix<int>& i, const Matrix<int>& j) const{ 
      return (*this)(i,j);
    }
    const Matrix<T> indexed(const CRSSparsity &sp) const{ return (*this)(sp); }
    
    /// set a non-zero
    void indexed_one_based_assignment(int k, const T & m){ at(k-1) = m;}
    void indexed_zero_based_assignment(int k, const T & m){ at(k) = m;}
    void indexed_assignment(const Slice &k, const Matrix<T>& m){ (*this)[k] = m;}
    void indexed_one_based_assignment(const Matrix<int> &k, const Matrix<T>& m){ (*this)[k-1] = m;}
    void indexed_zero_based_assignment(const Matrix<int> &k, const Matrix<T>& m){ (*this)[k] = m;}
    void indexed_assignment(const IndexList &k, const Matrix<T>& m){
      (*this)[k.getAll(size())] = m;
    }
    
    /// set a matrix element
    void indexed_one_based_assignment(int i, int j, const T & m){ elem(i-1,j-1) = m;}
    void indexed_zero_based_assignment(int i, int j, const T & m){ elem(i,j) = m;}
    void indexed_assignment(const Slice &i, const Slice &j, const Matrix<T>& m){ (*this)(i,j) = m; }
    void indexed_assignment(const IndexList &i, const IndexList &j, const Matrix<T>& m){
      (*this)(i.getAll(size1()),j.getAll(size2())) = m;
    }
    void indexed_assignment(const Slice &i, const Matrix<int>& j, const Matrix<T>& m){
      (*this)(i,j) = m;
    }
    void indexed_assignment( const Matrix<int>& i, const Slice &j, const Matrix<T>& m){
      (*this)(i,j) = m;
    }
    void indexed_assignment(const IndexList &i, const Matrix<int>& j, const Matrix<T>& m){
      (*this)(i.getAll(size1()),j) = m;
    }
    void indexed_assignment( const Matrix<int>& i, const IndexList &j, const Matrix<T>& m){
      (*this)(i,j.getAll(size2())) = m;
    } 
    void indexed_assignment( const Matrix<int>& i, const Matrix<int>& j, const Matrix<T>& m){
      (*this)(i,j) = m;
    } 
    void indexed_assignment(const CRSSparsity &sp,const Matrix<T>& m){
      (*this)(sp) = m;
    }
    //@}
    
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

  //@{
    /** \brief  Create nodes by their ID */
    static Matrix<T> binary(int op, const Matrix<T> &x, const Matrix<T> &y);
    static Matrix<T> unary(int op, const Matrix<T> &x);
    static Matrix<T> scalar_matrix(int op, const Matrix<T> &x, const Matrix<T> &y);
    static Matrix<T> matrix_scalar(int op, const Matrix<T> &x, const Matrix<T> &y);
    static Matrix<T> matrix_matrix(int op, const Matrix<T> &x, const Matrix<T> &y);
  //@}
    static void binary_no_alloc(void (*fcn)(unsigned char op, const T&, const T&, T&), unsigned char op, const Matrix<T> &x, const Matrix<T> &y, Matrix<T>& r, const std::vector<unsigned char>& mapping);

    /** \brief  Unary function */
#ifndef SWIG
    Matrix<T> unary_old(T (*fcn)(const T&)) const;
    Matrix<T> binary_old(T (*fcn)(const T&, const T&), const Matrix<T>& y) const;
        
    void unary_old(T (*fcn)(const T&), const Matrix<T>& x);
    void binary_old(T (*fcn)(const T&, const T&), const Matrix<T> &x, const Matrix<T> &y);
    void matrix_matrix_old(T (*fcn)(const T&, const T&), const Matrix<T>& x, const Matrix<T>& y);
    void matrix_scalar_old(T (*fcn)(const T&, const T&), const Matrix<T>& x, const T& y);
    void scalar_matrix_old(T (*fcn)(const T&, const T&), const T& x, const Matrix<T>& y);
#endif

    //@{
    /// Elementwise operations -- Octave/Python naming
    Matrix<T> __add__(const Matrix<T> &y) const;
    Matrix<T> __sub__(const Matrix<T> &y) const;
    Matrix<T> __mul__(const Matrix<T> &y) const;
    Matrix<T> __div__(const Matrix<T> &y) const;
    Matrix<T> __pow__(const Matrix<T> &y) const;
    Matrix<T> __constpow__(const Matrix<T> &y) const;
    Matrix<T> __mpower__(const Matrix<T> &y) const;
    Matrix<T> __mrdivide__  (const Matrix<T> &y) const;
    Matrix<T> __mldivide__   (const Matrix<T> &y) const;
    //@}
    
    /// Matrix product
    Matrix<T> mul(const Matrix<T> &y) const;

    /// Matrix product, no memory allocation: z += mul(x,y)
    static void mul_no_alloc(const Matrix<T> &x, const Matrix<T> &y_trans, Matrix<T>& z);

    /// Matrix product, no memory allocation: x += mul(z,trans(y))
    static void mul_no_alloc1(Matrix<T> &x, const Matrix<T> &y_trans, const Matrix<T>& z);

    /// Matrix product, no memory allocation: y += mul(trans(x),z)
    static void mul_no_alloc2(const Matrix<T> &x, Matrix<T> &y_trans, const Matrix<T>& z);
    
    /// Propagate sparsity using 0-1 logic through a matrix product, no memory allocation: z = mul(x,y)
    static void mul_sparsity(Matrix<T> &x, Matrix<T> &y_trans, Matrix<T>& z, bool fwd);
    
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
    Matrix<T> sign() const;
    Matrix<T> erfinv() const;
    Matrix<T> fmin(const Matrix<T>& y) const;
    Matrix<T> fmax(const Matrix<T>& y) const;
    Matrix<T> erf() const;
    Matrix<T> sinh() const;
    Matrix<T> cosh() const;
    Matrix<T> tanh() const;
    Matrix<T> log10() const;
    Matrix<T> printme(const Matrix<T>& y) const;
    //@}
    
    //@{
    /// Printing
#ifndef SWIG
    virtual void print(std::ostream &stream=std::cout) const; // print print description
    virtual void repr(std::ostream &stream=std::cout) const; // print representation
#endif
    static std::string className(); // name of the class
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
    
    /// Const access the sparsity - reference to data member
    const CRSSparsity& sparsity() const{ return sparsity_; }
    
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

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation 
    kl:    The number of subdiagonals in res 
    ku:    The number of superdiagonals in res 
    ldres: The leading dimension in res 
    res:   The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, T *res) const;

    // get the number if non-zeros for a given sparsity pattern
    int size(Sparsity sp) const;
    
    //@{
    /** \brief  create a sparse matrix with all zeros */
    static Matrix<T> sparse(int nrow, int ncol=1);
    static Matrix<T> sparse(const std::pair<int,int>& nm);
    //@}

    //@{
    /** \brief  create a dense matrix with all zeros */
    static Matrix<T> zeros(int nrow, int ncol=1);
    static Matrix<T> zeros(const std::pair<int,int>& nm);
    //@}

    //@{
    /** \brief  create a matrix with all ones */
    static Matrix<T> ones(int nrow, int ncol=1);
    static Matrix<T> ones(const std::pair<int,int>& nm);
    //@}

    //@{
    /** \brief  create a matrix with all inf */
    static Matrix<T> inf(int nrow=1, int ncol=1);
    static Matrix<T> inf(const std::pair<int,int>& nm);
    //@}
    
    //@{
    /** \brief  create a matrix with all nan */
    static Matrix<T> nan(int nrow=1, int ncol=1);
    static Matrix<T> nan(const std::pair<int,int>& nm);
    //@}

    //@{
    /** \brief  create a matrix by repeating an existing matrix */
    static Matrix<T> repmat(const Matrix<T>& x, int nrow, int ncol=1);
    static Matrix<T> repmat(const Matrix<T>& x, const std::pair<int,int>& nm);
    //@}

    /** \brief  create an n-by-n identity matrix */
    static Matrix<T> eye(int nrow);

  private:
    /// Sparsity of the matrix in a compressed row storage (CRS) format
    CRSSparsity sparsity_;
    
    /// Nonzero elements
    std::vector<T> data_;
};

} // namespace CasADi

// Typedefs/template initializations
namespace CasADi{
  typedef Matrix<int> IMatrix;
  typedef Matrix<double> DMatrix;
  typedef std::vector<Matrix<double> > DMatrixVector;
  typedef std::vector< std::vector<Matrix<double> > > DMatrixVectorVector;

  typedef DMatrix* DMatrixPtr;
  typedef std::vector<DMatrixPtr> DMatrixPtrV;
  typedef std::vector<DMatrixPtrV> DMatrixPtrVV;
} // namespace CasADi

#include "matrix_impl.hpp"

#endif // MATRIX_HPP

