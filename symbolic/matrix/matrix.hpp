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
#include "../stl_vector_tools.hpp"
#include "generic_matrix.hpp"
#include "generic_expression.hpp"

namespace CasADi{

  template<class T>
  struct NonZero {
    int k; // Non-zero index into matrix
    int i; // Row into matrix
    int j; // Col into matrix
    T el;  // Element
  };

  template<class T>
  class NonZeroIterator : public std::iterator< std::forward_iterator_tag, NonZero<T> > {
  public:
    NonZeroIterator(const Matrix<T> & m);

#ifndef SWIG
    NonZeroIterator<T>& operator++();
#endif // SWIG

    NonZero<T>& operator*();
    
    NonZeroIterator<T> begin();
    NonZeroIterator<T> end();
    bool operator==(const NonZeroIterator<T>& rhs);

  private:
    Matrix<T> m_;
    NonZero<T> nz;
  };

    
  //@{
  /** \brief Get typename */
  template <typename T> inline std::string matrixName() { return std::string("Matrix<") + typeid(T).name() + std::string(">");}
  template<> inline std::string matrixName<double>() { return "DMatrix"; }
  template<> inline std::string matrixName<int>() { return "IMatrix"; }
  //@}

  /** \brief General sparse matrix class
      General sparse matrix class that is designed with the idea that "everything is a matrix", that is, also scalars and vectors.\n
      This philosophy makes it easy to use and to interface in particularily with Python and Matlab/Octave.\n
  
      The syntax tries to stay as close as possible to the ublas syntax  when it comes to vector/matrix operations.\n

      Index starts with 0.\n
      Index vec happens as follows: (rr,cc) -> k = rr+cc*size1()\n
      Vectors are column vectors.\n
  
      The storage format is Compressed Column Storage (CCS), similar to that used for sparse matrices in Matlab, \n
      but unlike this format, we do allow for elements to be structurally non-zero but numerically zero.\n
  
      Matrix<T> is polymorphic with a std::vector<T> that contain all non-identical-zero elements.\n
      The sparsity can be accessed with Sparsity& sparsity()\n
  
      \author Joel Andersson 
      \date 2010        
  */
  template<class T>
  class Matrix : public GenericExpression<Matrix<T> >, public GenericMatrix<Matrix<T> >, public PrintableObject{
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
    
    /// Empty n-by-m matrix constructor
    Matrix(int nrow, int ncol);
    
    /// Dense n-by-m matrix filled with val constructor
    Matrix(int nrow, int ncol, const T& val);

    /// Sparse n-by-m matrix filled with given sparsity
    Matrix(int nrow, int ncol, const std::vector<int>& colind, const std::vector<int>& row, const std::vector<T>& d=std::vector<T>());

    /// Dense matrix constructor with data given as vector of vectors
    explicit Matrix(const std::vector< std::vector<T> >& m);
    
    //@{
    /// Sparse matrix with a given sparsity
    explicit Matrix(const Sparsity& sparsity, const T& val=0);
    //@}
    
    /// Sparse matrix with a given sparsity and non-zero elements.
    Matrix(const Sparsity& sparsity, const std::vector<T>& d);
    
    /** \brief Check if the dimensions and colind,row vectors are compatible.
     * \param complete  set to true to also check elementwise
     * throws an error as possible result
     */
    void sanityCheck(bool complete=false) const;
    
    /// This constructor enables implicit type conversion from a numeric type
    Matrix(double val);

    /// Construct from a vector
    /**
     * Thanks to implicit conversion, you can pretend that Matrix(const SXElement& x); exists.
     * Note: above remark applies only to C++, not python or octave interfaces
     */
    Matrix(const std::vector<T>& x);
    
    /// Construct dense matrix from a vector with the elements in column major ordering
    Matrix(const std::vector<T>& x, int nrow, int ncol);

    /// Convert to scalar type
    const T toScalar() const;
    
    /// Scalar type
    typedef T ScalarType;
    
#ifndef SWIG
    /// Base class
    typedef GenericMatrix<Matrix<T> > B;

    /// Expose base class functions
    using B::size;
    using B::sizeL;
    using B::sizeU;
    using B::numel;
    using B::size1;
    using B::size2;
    using B::shape;
    using B::empty;
    using B::scalar;
    using B::dense;
    using B::vector;
    using B::dimString;
    using B::sym;
    using B::zeros;
    using B::ones;
    using B::operator[];
    using B::operator();

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
    Matrix(const Matrix<A>& x) : sparsity_(x.sparsity()), data_(std::vector<T>(x.size())){
      copy(x.begin(),x.end(),begin());
    }

    /** \brief  Create an expression from an stl vector  */
    template<typename A>
    Matrix(const std::vector<A>& x) : sparsity_(Sparsity(x.size(),1,true)), data_(std::vector<T>(x.size())){
      copy(x.begin(),x.end(),begin());
    }

    /** \brief  Create a non-vector expression from an stl vector */
    template<typename A>
    Matrix(const std::vector<A>& x,  int nrow, int ncol) : sparsity_(Sparsity(nrow,ncol,true)), data_(std::vector<T>(x.size())){
      if(x.size() != nrow*ncol) throw CasadiException("Matrix::Matrix(const std::vector<T>& x,  int n, int m): dimension mismatch");
      copy(x.begin(),x.end(),begin());
    }
    
    /** \brief  ublas vector */
#ifdef HAVE_UBLAS
    template<typename T, typename A>
    explicit Matrix<T>(const ublas::vector<A> &x) : sparsity_(Sparsity(x.size(),1,true)), data_(std::vector<T>(x.size())){
      copy(x.begin(),x.end(),begin());
    }

    template<typename T, typename A>
    explicit Matrix<T>(const ublas::matrix<A> &x) : sparsity_(Sparsity(x.size1(),x.size2(),true)), data_(std::vector<T>(numel())){
      copy(x.begin(),x.end(),begin());
      return ret;
    }
#endif // HAVE_UBLAS
#endif // SWIG


#ifndef SWIG
    /// Get a non-zero element
    inline const T& at(int k) const{
      return const_cast<Matrix<T>*>(this)->at(k);
    }
    
    /// Access a non-zero element
    inline T& at(int k){
      try{
        if (k<0) k+=size(); 
        return data().at(k);
      } catch(std::out_of_range& ex){
        std::stringstream ss;
        ss << "Out of range error in Matrix<>::at: " << k << " not in range [0," << size() << ")";
        throw CasadiException(ss.str());
      }
    }
#else // SWIG
    /// Access a non-zero element
    T at(int k){
      try{
        if (k<0) k+=size(); 
        return data().at(k);
      } catch(std::out_of_range& ex){
        std::stringstream ss;
        ss << "Out of range error in Matrix<>::at: " << k << " not in range [0," << size() << ")";
        throw CasadiException(ss.str());
      }
    }
#endif // SWIG
    
#ifndef SWIG
    /// get an element
    const T& elem(int rr, int cc=0) const;
    
    /// get a reference to an element
    T& elem(int rr, int cc=0);
#else // SWIG
    /// Access a non-zero element
    T elem(int rr, int cc=0) { return elem(rr,cc);}
#endif // SWIG

    /// get an element, do not allocate
    const T getElement(int rr, int cc=0) const{ return elem(rr,cc);}
    
    /// Returns true if the matrix has a non-zero at location rr,cc
    bool hasNZ(int rr, int cc) const { return sparsity().hasNZ(rr,cc); }

    /// Returns the truth value of a Matrix
    bool __nonzero__() const;

    //@{
    /// Get a submatrix
    const Matrix<T> sub(int rr, int cc) const;
    const Matrix<T> sub(const std::vector<int>& rr, int cc) const{ return sub(rr,std::vector<int>(1,cc));}
    const Matrix<T> sub(int rr, const std::vector<int>& cc) const{ return sub(std::vector<int>(1,rr),cc);}
    const Matrix<T> sub(const std::vector<int>& rr, const std::vector<int>& cc) const;
    const Matrix<T> sub(const Slice& rr, const std::vector<int>& cc) const { return sub(rr.getAll(size1()),cc);}
    const Matrix<T> sub(const std::vector<int>& rr, const Slice& cc) const { return sub(rr,cc.getAll(size2()));}
    const Matrix<T> sub(const Slice& rr, const Slice& cc) const{ return sub(rr.getAll(size1()),cc.getAll(size2()));}
    const Matrix<T> sub(const Slice& rr, int cc) const{ return sub(rr.getAll(size1()),std::vector<int>(1,cc));}
    const Matrix<T> sub(int rr, const Slice& cc) const{ return sub(std::vector<int>(1,rr),cc.getAll(size2()));}
    const Matrix<T> sub(const Matrix<int>& rr, const std::vector<int>& cc) const;
    const Matrix<T> sub(const std::vector<int>& rr, const Matrix<int>& cc) const;
    const Matrix<T> sub(const Matrix<int>& rr, const Slice& cc) const {return sub(rr,cc.getAll(size2()));}
    const Matrix<T> sub(const Slice& rr, const Matrix<int>& cc) const {return sub(rr.getAll(size1()),cc);}
    const Matrix<T> sub(const Matrix<int>& rr, const Matrix<int>& cc) const;
    const Matrix<T> sub(const Sparsity& sp, int dummy = 0) const;
    //@}

    //@{
    /// Set a submatrix
    void setSub(const Matrix<T>& m, int rr, int cc);
    void setSub(const Matrix<T>& m, const std::vector<int>& rr, int cc){ setSub(m,rr,std::vector<int>(1,cc));}
    void setSub(const Matrix<T>& m, int rr, const std::vector<int>& cc){ setSub(m,std::vector<int>(1,rr),cc);}
    void setSub(const Matrix<T>& m, const std::vector<int>& rr, const std::vector<int>& cc);
    void setSub(const Matrix<T>& m, const Slice& rr, const std::vector<int>& cc){ setSub(m,rr.getAll(size1()),cc);}
    void setSub(const Matrix<T>& m, const std::vector<int>& rr, const Slice& cc){ setSub(m,rr,cc.getAll(size2()));}
    void setSub(const Matrix<T>& m, const Slice& rr, const Slice& cc){ setSub(m,rr.getAll(size1()),cc.getAll(size2()));}
    void setSub(const Matrix<T>& m, const Matrix<int>& rr, const std::vector<int>& cc);
    void setSub(const Matrix<T>& m, const std::vector<int>& rr, const Matrix<int>& cc);
    void setSub(const Matrix<T>& m, const Matrix<int>& rr, const Slice& cc) {return setSub(m,rr,cc.getAll(size2()));}
    void setSub(const Matrix<T>& m, const Slice& rr, const Matrix<int>& cc) {return setSub(m,rr.getAll(size1()),cc);}
    void setSub(const Matrix<T>& m, const Matrix<int>& rr, const Matrix<int>& cc);
    void setSub(const Matrix<T>& m, const Sparsity& sp, int dummy);
    //@}

    //@{
    /// Add a submatrix to an existing matrix (TODO: remove memory allocation)
    template<typename RR, typename CC>
    void addSub(const Matrix<T>& m, RR rr, CC cc){ setSub(m+sub(rr,cc),rr,cc);}
    //@}

    //@{
    /// Retrieve a submatrix (TODO: remove memory allocation)
    template<typename RR, typename CC>
    void getSub(Matrix<T>& m, RR rr, CC cc){ m = sub(rr,cc);}
    //@}

    //@{
    /// Get a set of nonzeros
    const Matrix<T> getNZ(int k) const{ return at(k);}
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

    /// Append a matrix vertically (NOTE: only efficient if vector)
    void append(const Matrix<T>& y);

    /// Append a matrix horizontally
    void appendColumns(const Matrix<T>& y);

    //@{
    /// Indexing for interfaced languages
    /// get a non-zero
    const Matrix<T> indexed_one_based(int k) const{ return this->operator[](k-1);}
    const Matrix<T> indexed_zero_based(int k) const{ return this->operator[](k);}
    const Matrix<T> indexed_one_based(const Matrix<int>& k) const{ return this->operator[](k-1);}
    const Matrix<T> indexed_zero_based(const Matrix<int>& k) const{ return this->operator[](k);}
    const Matrix<T> indexed(const Slice &k) const{ return this->operator[](k);}
    const Matrix<T> indexed(const IndexList &k) const{
      return (*this)[k.getAll(size())];
    }
    
    /// get a matrix element
    const Matrix<T> indexed_one_based(int rr, int cc) const{ return (*this)(rr-1,cc-1);}
    const Matrix<T> indexed_zero_based(int rr, int cc) const{ return (*this)(rr,cc);}
    const Matrix<T> indexed(const Slice &rr, const Slice &cc) const{ return (*this)(rr,cc); }
    const Matrix<T> indexed(const IndexList &rr, const IndexList &cc) const{ 
      return (*this)(rr.getAll(size1()),cc.getAll(size2()));
    }
    const Matrix<T> indexed(const Slice &rr, const Matrix<int>& cc) const{ return (*this)(rr,cc); }
    const Matrix<T> indexed(const Matrix<int>& rr, const IndexList &cc) const{ 
      return (*this)(rr,cc.getAll(size2()));
    }
    const Matrix<T> indexed(const Matrix<int>& rr, const Slice &cc) const{ return (*this)(rr,cc); }
    const Matrix<T> indexed(const IndexList& rr, const Matrix<int> &cc) const{ 
      return (*this)(rr.getAll(size1()),cc);
    }
    const Matrix<T> indexed(const Matrix<int>& rr, const Matrix<int>& cc) const{ 
      return (*this)(rr,cc);
    }
    const Matrix<T> indexed(const Sparsity &sp) const{ return (*this)(sp); }
    
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
    void indexed_one_based_assignment(int rr, int cc, const T & m){ elem(rr-1,cc-1) = m;}
    void indexed_zero_based_assignment(int rr, int cc, const T & m){ elem(rr,cc) = m;}
    void indexed_assignment(const Slice &rr, const Slice &cc, const Matrix<T>& m){ (*this)(rr,cc) = m; }
    void indexed_assignment(const IndexList &rr, const IndexList &cc, const Matrix<T>& m){
      (*this)(rr.getAll(size1()),cc.getAll(size2())) = m;
    }
    void indexed_assignment(const Slice &rr, const Matrix<int>& cc, const Matrix<T>& m){
      (*this)(rr,cc) = m;
    }
    void indexed_assignment( const Matrix<int>& rr, const Slice &cc, const Matrix<T>& m){
      (*this)(rr,cc) = m;
    }
    void indexed_assignment(const Matrix<int> &rr, const IndexList& cc, const Matrix<T>& m){
      (*this)(rr,cc.getAll(size2())) = m;
    }
    void indexed_assignment( const IndexList& rr, const Matrix<int> &cc, const Matrix<T>& m){
      (*this)(rr.getAll(size1()),cc) = m;
    } 
    void indexed_assignment( const Matrix<int>& rr, const Matrix<int>& cc, const Matrix<T>& m){
      (*this)(rr,cc) = m;
    } 
    void indexed_assignment(const Sparsity &sp,const Matrix<T>& m){
      // (*this)(sp) = m;   // VC2010 compiler errors
          SubMatrix<Matrix<T>,Sparsity,int> temp(*this,sp,0);
          temp = m;
    }
    //@}
    
    /// Set all elements to zero
    void setZero();
    
    /// Set all elements to a value
    void setAll(const T& val);
    
    /** \brief  Make the matrix an dense nrow-by-ncol matrix */
    void makeDense(int nrow, int ncol, const T& val);

    /** \brief  Make the matrix an empty nrow-by-ncol matrix */
    void makeEmpty(int nrow, int ncol);

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
  
    //@{
    /// Elementwise operations -- Octave/Python naming
    Matrix<T> __add__(const Matrix<T> &y) const;
    Matrix<T> __sub__(const Matrix<T> &y) const;
    Matrix<T> __mul__(const Matrix<T> &y) const;
    Matrix<T> __div__(const Matrix<T> &y) const;
    Matrix<T> __lt__(const Matrix<T> &y) const;
    Matrix<T> __le__(const Matrix<T> &y) const;
    Matrix<T> __eq__(const Matrix<T> &y) const;
    Matrix<T> __ne__(const Matrix<T> &y) const;
    Matrix<T> __truediv__(const Matrix<T> &y) const {return __div__(y);};
    Matrix<T> __pow__(const Matrix<T> &y) const;
    Matrix<T> __constpow__(const Matrix<T> &y) const;
    Matrix<T> __mpower__(const Matrix<T> &y) const;
    Matrix<T> __mrdivide__  (const Matrix<T> &y) const;
    //@}
    
    /// Matrix-matrix product
    Matrix<T> mul_full(const Matrix<T> &y, const Sparsity & sp_z=Sparsity()) const;

    /// Matrix-matrix product
    Matrix<T> mul(const Matrix<T> &y, const Sparsity & sp_z=Sparsity()) const;
    
    /// Matrix-matrix product, no memory allocation: z += mul(x,y)
    static void mul_no_alloc_nn(const Matrix<T> &x, const Matrix<T>& y, Matrix<T>& z);
    
    /// Matrix-matrix product, no memory allocation: z += mul(trans(x),y)
    static void mul_no_alloc_tn(const Matrix<T> &trans_x, const Matrix<T> &y, Matrix<T>& z);

    /// Matrix-matrix product, no memory allocation: z += mul(x,trans(y))
    static void mul_no_alloc_nt(const Matrix<T>& x, const Matrix<T> &trans_y, Matrix<T>& z);
  
    /// Matrix-vector product, no memory allocation: z += mul(trans(x),y)
    static void mul_no_alloc_tn(const Matrix<T>& trans_x, const std::vector<T> &y, std::vector<T>& z);

    /// vector-matrix product, no memory allocation: z += mul(x,y)
    static void mul_no_alloc_nn(const Matrix<T>& x, const std::vector<T> &y, std::vector<T>& z);
  
    /// Propagate sparsity using 0-1 logic through a matrix product, no memory allocation: z = mul(trans(x),y)
    template<bool Fwd>
    static void mul_sparsity(Matrix<T> &x_trans, Matrix<T> &y, Matrix<T>& z);
  
    /// Calculates inner_prod(x,mul(A,x)) without memory allocation
    static T quad_form(const Matrix<T>& A, const std::vector<T>& x);
  
    /// Matrix transpose
    Matrix<T> trans() const;
    
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
    Matrix<T> __copysign__(const Matrix<T>& y) const;
    Matrix<T> erfinv() const;
    Matrix<T> fmin(const Matrix<T>& y) const;
    Matrix<T> fmax(const Matrix<T>& y) const;
    Matrix<T> erf() const;
    Matrix<T> sinh() const;
    Matrix<T> cosh() const;
    Matrix<T> tanh() const;
    Matrix<T> arcsinh() const;
    Matrix<T> arccosh() const;
    Matrix<T> arctanh() const;
    Matrix<T> arctan2(const Matrix<T>& y) const;
    Matrix<T> log10() const;
    Matrix<T> printme(const Matrix<T>& y) const;
    Matrix<T> logic_not() const;
    Matrix<T> logic_and(const Matrix<T>& y) const;
    Matrix<T> logic_or(const Matrix<T>& y) const;
    Matrix<T> if_else_zero(const Matrix<T>& y) const;
    //@}
    
    //@{
    /// Printing
#ifndef SWIG
    virtual void print(std::ostream &stream=std::cout) const; // print print description
    virtual void repr(std::ostream &stream=std::cout) const; // print representation
#endif
    static std::string className(); // name of the class
    void printScalar(std::ostream &stream=std::cout) const; // print scalar
    void printVector(std::ostream &stream=std::cout) const; // print vector-style
    void printDense(std::ostream &stream=std::cout) const; // Print dense matrix-stype
    void printSparse(std::ostream &stream=std::cout) const; // print sparse matrix style
    //@}
  
    // Get the sparsity pattern
    const std::vector<int>& row() const;
    const std::vector<int>& colind() const;
    int row(int el) const;
    int colind(int col) const;
    void clear();
    void resize(int nrow, int ncol);
    void reserve(int nnz);
    void reserve(int nnz, int ncol);
    
    /** \brief Erase a submatrix
        Erase rows and/or columns of a matrix */
    void erase(const std::vector<int>& rr, const std::vector<int>& cc);
    
    /** \brief Remove cols or rows
        Rremove/delete rows and/or columns of a matrix */
    void remove(const std::vector<int>& rr, const std::vector<int>& cc);
    
    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros */
    void enlarge(int nrow, int ncol, const std::vector<int>& rr, const std::vector<int>& cc);
    
    /// Access the non-zero elements
    std::vector<T>& data();
    
    /// Const access the non-zero elements
    const std::vector<T>& data() const;
    
    /// Get a pointer to the data
    T* ptr(){ return empty() ? static_cast<T*>(0) : &front();}
    
    /// Get a const pointer to the data
    const T* ptr() const{ return empty() ? static_cast<const T*>(0) : &front();}
        
    /// Const access the sparsity - reference to data member
    const Sparsity& sparsity() const{ return sparsity_; }
    
    /// Access the sparsity, make a copy if there are multiple references to it
    Sparsity& sparsityRef();
    
    /** \brief  Set the non-zero elements, scalar */
    void set(T val, SparsityType sp=SPARSE);
    
    /** \brief  Get the non-zero elements, scalar */
    void get(T& val, SparsityType sp=SPARSE) const;

    /** \brief  Set the non-zero elements, vector */
    void set(const std::vector<T>& val, SparsityType sp=SPARSE);

    /** \brief  Get the non-zero elements, vector */
    void get(std::vector<T>& val, SparsityType sp=SPARSE) const;

    /** \brief  Set the non-zero elements, Matrix */
    void set(const Matrix<T>& val, SparsityType sp=SPARSE);

    /** \brief  Get the non-zero elements, Matrix */
    void get(Matrix<T>& val, SparsityType sp=SPARSE) const;

#ifdef SWIG
    %rename(get) getStridedArray;
    %rename(set) setArray;
#endif

    /** \brief  Get the non-zero elements, array */
    void getArray(T* val, int len, SparsityType sp=SPARSE) const;

    /** \brief  Set the non-zero elements, array */
    void setArray(const T* val, int len, SparsityType sp=SPARSE);

    /** \brief  Get the non-zero elements, array, sparse and correct length */
    void getArray(T* val) const;

    /** \brief  Set the non-zero elements, array, sparse and correct length */
    void setArray(const T* val);
    
    /** \brief  Get the non-zero elements, strided array */
    void getStridedArray(T* val, int len, int stride1, int stride2, SparsityType sp=SPARSE) const;
    
#ifndef SWIG
    /** \brief  Legacy - use getArray instead */
    void get(T* val, SparsityType sp=SPARSE) const;

    /** \brief  Legacy - use setArray instead */
    void set(const T* val, SparsityType sp=SPARSE);

    /** Bitwise set, reinterpreting the data as a bvec_t array */
    void setZeroBV();

    /** Bitwise set, reinterpreting the data as a bvec_t array */
    void setBV(const Matrix<T>& val);

    /** Bitwise set, reinterpreting the data as a bvec_t array */
    void getBV(Matrix<T>& val) const{ val.setBV(*this);}

    /** Bitwise or, reinterpreting the data as a bvec_t array */
    void borBV(const Matrix<T>& val);

    /** \brief Bitwise get the non-zero elements, array */
    void getArrayBV(bvec_t* val, int len) const;

    /** \brief Bitwise set the non-zero elements, array */
    void setArrayBV(const bvec_t* val, int len);

    /** \brief Bitwise or the non-zero elements, array */
    void borArrayBV(const bvec_t* val, int len);

#endif

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation 
        kl:    The number of subdiagonals in res 
        ku:    The number of superdiagonals in res 
        ldres: The leading dimension in res 
        res:   The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, T *res) const;
    
    //@{
    /** \brief  create a sparse matrix with all zeros */
    static Matrix<T> sparse(int nrow=1, int ncol=1);
    static Matrix<T> sparse(const std::pair<int,int>& rc){ return sparse(rc.first,rc.second);}
    //@}
    
    /* \brief Construct a sparse matrix from triplet form
     * Default matrix size is max(col) x max(row)
     */
    //@{
    static Matrix<T> sparse(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d);
    static Matrix<T> sparse(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d, int nrow, int ncol);
    static Matrix<T> sparse(const std::vector<int>& row, const std::vector<int>& col, const std::vector<T>& d, const std::pair<int,int>& rc);
    //@}
    
    //@{
    /** \brief  create a matrix with all inf */
    static Matrix<T> inf(const Sparsity& sp);
    static Matrix<T> inf(int nrow=1, int ncol=1);
    static Matrix<T> inf(const std::pair<int,int>& rc);
    //@}
    
    //@{
    /** \brief  create a matrix with all nan */
    static Matrix<T> nan(const Sparsity& sp);
    static Matrix<T> nan(int nrow=1, int ncol=1);
    static Matrix<T> nan(const std::pair<int,int>& rc);
    //@}

    //@{
    /** \brief  create a matrix by repeating an existing matrix */
    static Matrix<T> repmat(const T& x, const Sparsity& sp);
    static Matrix<T> repmat(const Matrix<T>& x, const Sparsity& sp);
    static Matrix<T> repmat(const Matrix<T>& x, int nrow, int ncol=1);
    static Matrix<T> repmat(const Matrix<T>& x, const std::pair<int,int>& rc);
    //@}

    /** \brief  create an n-by-n identity matrix */
    static Matrix<T> eye(int ncol);

    /** \brief  The following function is used to ensure similarity to MX, which is reference counted */
    bool isNull() const{ return false;}
    
    // @{
    /// Set the 'precision, width & scientific' used in printing and serializing to streams
    static void setPrecision(int precision) { stream_precision_ = precision; }
    static void setWidth(int width) { stream_width_ = width; }
    static void setScientific(bool scientific) { stream_scientific_ = scientific; }
    // @}
    
  private:
    /// Sparsity of the matrix in a compressed column storage (CCS) format
    Sparsity sparsity_;
    
    /// Nonzero elements
    std::vector<T> data_;
    
    /// Precision used in streams
    static int stream_precision_;
    static int stream_width_;
    static bool stream_scientific_;
    
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

