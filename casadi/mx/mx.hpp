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

#ifndef MX_HPP
#define MX_HPP

#include "../shared_object.hpp"
#include "../matrix/matrix.hpp"
#include "../sx/sx.hpp"
#include <vector>
namespace CasADi{
  
/** \brief  Forward declaration */
class MXNode;


/** \brief MX - Matrix expression
  The MX class is used to build up trees made up from MXNodes. It is a more general graph representation than the scalar expression,
  SX, and much less efficient for small objects. On the other hand, the class allows much more general operations than does SX,
  in particular matrix valued operations and calls to arbitrary differentiable functions.
  
  The MX class is designed to have identical syntax with the Matrix<> template class, and uses Matrix<double> as its internal 
  representation of the values at a node. By keeping the syntaxes identical, it is possible to switch from one class to the other, 
  as well as inlining MX functions to SX functions.
  
  Note that an operation is always "lazy", making a matrix multiplication will create a matrix multiplication node, not perform
  the actual multiplication.

  \author Joel Andersson 
  \date 2010-2011
*/
class MX : public SharedObject{
  public:
  
    /** \brief  Default constructor */
    MX();

    /** \brief  Construct a symbolic matrix (matrix variable) */
    explicit MX(const std::string& name, int n=1, int m=1);

    /** \brief  Construct a symbolic matrix (matrix variable) */
    explicit MX(const std::string& name, const std::pair<int,int> &nm);

    /** \brief  Construct a symbolic matrix (matrix variable) */
    explicit MX(const std::string& name, const CRSSparsity & sp);

    /** \brief  Construct MX with a given sparsity */
    explicit MX(const CRSSparsity& sp, const MX& val);
    
    /** \brief  Create scalar constant (also implicit type conversion) */
    MX(double x);

    /** \brief  Copy constructor */
    MX(const MX& x);

    /** \brief  Create vector constant (also implicit type conversion) */
    MX(const std::vector<double> &x);
    
    /** \brief  Create sparse matrix constant (also implicit type conversion) */
    MX(const Matrix<double> &x);

    /** \brief  Matrix with all zeros */
    MX(int nrow, int ncol);
    
    /** \brief  Dense matrix filled with value val */
    MX(int nrow, int ncol, const MX& val);
    
    /** \brief  Destructor */
    virtual ~MX();

#ifndef SWIG
   /** \brief  Create from node */
    static MX create(MXNode* node);

    /** \brief  Get vector nonzero or slice of nonzeros */
    template<typename K>
    const MX operator[](const K& k) const{ return getNZ(k); }

    /** \brief  Access vector nonzero or slice of nonzeros */
    template<typename K>
    NonZeros<MX,K> operator[](const K& k){ return NonZeros<MX,K>(*this,k); }

    /** \brief  Get vector element or slice */
    template<typename I>
    const MX operator()(const I& i) const{ return getSub(i,0);}
    
    /** \brief  Get Matrix element or slice */
    template<typename I, typename J>
    const MX operator()(const I& i, const J& j) const{ return getSub(i,j); }

    /** \brief  Access vector element or slice */
    template<typename I>
    SubMatrix<MX,I,int> operator()(const I& i){ return SubMatrix<MX,I,int>(*this,i,0); }
    
    /** \brief  Access Matrix element or slice */
    template<typename I, typename J>
    SubMatrix<MX,I,J> operator()(const I& i, const J& j){ return SubMatrix<MX,I,J>(*this,i,j); }

    /// Get a non-zero element, with bounds checking
    const MX at(int k) const;

    /// Access a non-zero element, with bounds checking
    NonZeros<MX,int> at(int k);
    
#endif // SWIG
    
    //@{
    /// Indexing for interfaced languages
    
    /// get a non-zero
    const MX indexed_one_based(int k) const{ return at(k-1);}
    const MX indexed_zero_based(int k) const{ return at(k);}
    const MX indexed(const IndexList &k) const{
      return (*this)[k.getAll(size())];
    }
    const MX indexed(const Slice &k) const{ 
      return (*this)[k.getAll(size())];
    }
    
    /// get a matrix element
    const MX indexed_one_based(int i, int j) const{ return (*this)(i-1,j-1);}
    const MX indexed_zero_based(int i, int j) const{ return (*this)(i,j);}
    const MX indexed(const IndexList &i, const IndexList &j) const{ 
      return (*this)(i.getAll(size1()),j.getAll(size2()));
    }
    const MX indexed(const Slice &i, const Slice &j) const{ 
      return (*this)(i.getAll(size1()),j.getAll(size2()));
    }
    
    /// set a non-zero
    void indexed_one_based_assignment(int k, const MX &m){ at(k-1) = m(0,0);}
    void indexed_zero_based_assignment(int k, const MX &m){ at(k) = m(0,0);}
    void indexed_assignment(const IndexList &k, const MX &m){
      (*this)[k.getAll(size())] = m;
    }
    void indexed_assignment(const Slice &k, const MX &m){
      (*this)[k.getAll(size())] = m;
    }
    
    /// set a matrix element
    void indexed_one_based_assignment(int i, int j, const MX &m){ (*this)(i-1,j-1) = m;}
    void indexed_zero_based_assignment(int i, int j, const MX &m){ (*this)(i,j) = m;}
    void indexed_assignment(const IndexList &i, const IndexList &j, const MX &m){
      setSub(i.getAll(size1()),j.getAll(size2()),m);
    }
    
    void indexed_assignment(const Slice &i, const Slice &j, const MX &m){
      (*this)(i.getAll(size1()),j.getAll(size2())) = m;
    }
    //@}
    
    /// Scalar type
    typedef MX ScalarType;

    /** \brief  Get the number of (structural) non-zero elements */
    int size() const;

    /** \brief  Get the number of elements */
    int numel() const;

    /** \brief get the first dimension (i.e. n for a n-by-m matrix) */
    int size1() const;
    
    /** \brief get the first dimension (i.e. m for a n-by-m matrix) */
    int size2() const;

    /** \brief Get the sparsity pattern */
    const CRSSparsity& sparsity() const;

    /// Access the sparsity, make a copy if there are multiple references to it
    CRSSparsity& sparsityRef();

    /** \brief Erase a submatrix */
    void erase(const std::vector<int>& ii, const std::vector<int>& jj);

    /** \brief Enlarge matrix
    Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros */
    void enlarge(int nrow, int ncol, const std::vector<int>& ii, const std::vector<int>& jj);

  #ifndef SWIG  
  //@{
  /** \brief  Operators that changes the object */
    MX& operator+=(const MX &y);
    MX& operator-=(const MX &y);
    MX& operator*=(const MX &y);
    MX& operator/=(const MX &y);
  //@}


  //@{
  /** \brief  Operators */
    friend MX operator+(const MX &x, const MX &y);
    friend MX operator-(const MX &x, const MX &y);
    friend MX operator*(const MX &x, const MX &y);
    friend MX operator/(const MX &x, const MX &y);

  //@}
  #endif // SWIG
  
  MX operator-() const;

  /** \brief  Check if the matrix expression is empty */
  bool empty() const;
  
  /** \brief  Check if the matrix expression is dense */
  bool dense() const;
  
  //@{
  /** \brief  Access a member of the node */
  MXNode* operator->();

  /** \brief  Const access a member of the node */
  const MXNode* operator->() const;
  //@}

  //@{
    /** \brief  Create nodes by their ID */
    static MX binary(int op, const MX &x, const MX &y);
    static MX unary(int op, const MX &x);
    static MX scalar_matrix(int op, const MX &x, const MX &y);
    static MX matrix_scalar(int op, const MX &x, const MX &y);
    static MX matrix_matrix(int op, const MX &x, const MX &y);
  //@}

  /** \brief  Matrix of all zeros */  
  static MX zeros(int nrow, int ncol=1);
  static MX zeros(const CRSSparsity& sparsity);
  
  /** \brief  Matrix of all ones */  
  static MX ones(int nrow, int ncol=1);
  
  /** \brief  Identity matrix */  
  static MX eye(int nrow);

  /** \brief  Get the jacobian of an function evaluation with respect to the iind-th argument */
  MX jac(int iind=0);

  const MX getSub(int i, int j) const;
  const MX getSub(int i, const std::vector<int>& j) const;
  const MX getSub(const std::vector<int>& i, int j) const;
  const MX getSub(const std::vector<int>& i, const std::vector<int>& j) const;
  
  void setSub(int i, int j, const MX& el);
  void setSub(int i, const std::vector<int>& j, const MX& el);
  void setSub(const std::vector<int>& i, int j, const MX& el);
  void setSub(const std::vector<int>& i, const std::vector<int>& j, const MX& el);
  
  MX getNZ(int k) const;
  MX getNZ(const std::vector<int>& k) const;
  void setNZ(int k, const MX& el);
  void setNZ(const std::vector<int>& k, const MX& el);
  
  /// Numeric evaluation
  Matrix<double> eval(const std::vector<Matrix<double> >& x);
  
  /// Symbolic evaluation (matrix graph)
  Matrix<SX> eval(const std::vector<Matrix<SX> >& x);
  
  /// Symbolic evaluation (matrix graph)
  MX eval(const std::vector<MX>& x);

  /** \brief Get string representation of dimensions.
  The representation is (nrow x ncol = numel | size)
  */
  std::string dimString() const;
  
  // all binary operations
  MX __add__(const MX& b) const;
  MX __radd__(const MX& b) const;
  MX __sub__(const MX& b) const;
  MX __rsub__(const MX& b) const;
  MX __mul__(const MX& b) const;
  MX __rmul__(const MX& b) const;
  MX __div__(const MX& b) const;
  MX __rdiv__(const MX& b) const;
  MX __pow__(const MX& b) const;
  MX __rpow__(const MX& b) const;
  MX __mrdivide__  (const MX& b) const;
  MX __rmrdivide__ (const MX& b) const;
  MX __ldivide__   (const MX& b) const;
  MX __rmldivide__ (const MX& b) const;
  MX __mpower__(const MX& b) const;
  MX __rmpower__(const MX& b) const;
  MX prod(const MX& y) const;
  MX rprod(const MX& y) const;
  MX inner_prod(const MX& y) const;
  MX outer_prod(const MX& y) const;
  MX constpow(const MX& y) const;
  MX fmin(const MX& y) const;
  MX fmax(const MX& y) const;

  // all unary operations
  MX exp() const;
  MX log() const;
  MX sqrt() const;
  MX sin() const;
  MX cos() const;
  MX tan() const;
  MX arcsin() const;
  MX arccos() const;
  MX arctan() const;
  MX floor() const;
  MX ceil() const;
  MX erf() const;
  MX sinh() const;
  MX cosh() const;
  MX tanh() const;

  /** \brief  Returns the IMatrix that represents the mapping of a Mapping node
  *
  */
  const Matrix<int>& mapping();

  
};

//@{
/// Some typedefs
typedef DMatrix* DMatrixPtr;
typedef std::vector<DMatrixPtr> DMatrixPtrV;
typedef std::vector<DMatrixPtrV> DMatrixPtrVV;

typedef SXMatrix* SXMatrixPtr;
typedef std::vector<SXMatrixPtr> SXMatrixPtrV;
typedef std::vector<SXMatrixPtrV> SXMatrixPtrVV;

typedef MX* MXPtr;
typedef std::vector<MXPtr> MXPtrV;
typedef std::vector<MXPtrV> MXPtrVV;
//@}

#ifdef SWIG 
%extend MX {
  std::string __repr__() { return $self->getRepresentation();}
}
#endif // SWIG
} // namespace CasADi

#ifdef SWIG 
// Template instantiations
%template(MXVector) std::vector<CasADi::MX>;
#endif // SWIG

#endif // MX_HPP
