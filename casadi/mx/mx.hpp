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
class MX;

}

#ifndef SWIG  
namespace std{
//@{
/** \brief  Functions with c equivalents: The implementation and syntax mirrors the standard c functions in math.h */
#define MX CasADi::MX
MX sqrt(const MX &x);
MX sin(const MX &x);
MX cos(const MX &x);
MX tan(const MX &x);
MX atan(const MX &x);
MX asin(const MX &x);
MX acos(const MX &x);
MX exp(const MX &x);
MX log(const MX &x);
MX pow(const MX &x, const MX &n);
MX constpow(const MX &x, const MX &n);
MX abs(const MX &x);
MX fabs(const MX &x); // same as abs
MX floor(const MX &x);
MX ceil(const MX &x);
MX erf(const MX &x);
MX fmin(const MX &a, const MX &b);
MX fmax(const MX &a, const MX &b);
#undef MX
//@}
} // namespace std
#endif // SWIG  

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
  static MX zeros(int nrow, int ncol);
  static MX zeros(const CRSSparsity& sparsity);
  
  /** \brief  Matrix of all ones */  
  static MX ones(int nrow, int ncol);
  
  /** \brief  Identity matrix */  
  static MX eye(int n);

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
  
  // all binary operations with a particular right argument
  #define binops(t) \
  MX __add__(t b){    return *this + b;} \
  MX __radd__(t b){   return b + *this;} \
  MX __sub__(t b){    return *this - b;} \
  MX __rsub__(t b){   return b - *this;} \
  MX __mul__(t b){    return *this * b;} \
  MX __rmul__(t b){   return b * *this;} \
  MX __div__(t b){    return *this / b;} \
  MX __rdiv__(t b){   return b / *this;} \
  MX __pow__(t b) const {    return std::pow(*this,b);} \
  MX __rpow__(t b) const {   return std::pow(b,*this);} \
  MX __mrdivide__  (t b) const { if (MX(b).numel()==1) return *this/b; throw CasadiException("mrdivide: Not implemented");} \
  MX __rmrdivide__ (t b) const { if ((*this).numel()==1) return b/(*this); throw CasadiException("rmrdivide: Not implemented");} \
  MX __ldivide__   (t b) const { if (MX(b).numel()==1) return *this/b; throw CasadiException("mldivide: Not implemented");} \
  MX __rmldivide__ (t b) const { if ((*this).numel()==1) return b/(*this); throw CasadiException("rmldivide: Not implemented");} \
  MX __mpower__(t b) const  {   return std::pow(*this,b); throw CasadiException("mpower: Not implemented");} \
  MX __rmpower__(t b) const {   return std::pow(b,*this); throw CasadiException("rmpower: Not implemented");}
  
  // Binary operations with all right hand sides
  binops(const MX&)
  binops(double)
  binops(const Matrix<double> &)
  #undef binops


  /** \brief  Returns the IMatrix that represents the mapping of a Mapping node
  *
  */
  const Matrix<int>& mapping();

};

//@{
/// Some typedefs
typedef double* Dptr;
typedef std::vector<double*> VDptr;
typedef std::vector<std::vector<double* > > VVDptr;
//@}


} // namespace CasADi

#ifdef SWIG 

namespace CasADi {

%extend MX {
std::string __repr__() { return $self->getRepresentation();
}

 // all binary operations with a particular right argument
  #define binops(t) \
  MX prod(t b){       return prod(*$self,b);} \
  MX rprod(t b){       return prod(b,*$self);} \
  MX inner_prod(t b){ return inner_prod(*$self,b);} \
  MX outer_prod(t b){ return outer_prod(*$self,b);} \
  MX constpow(t b){    return std::constpow(*$self,b);} \
  MX fmin(t b){       return std::fmin(*$self,b);} \
  MX fmax(t b){       return std::fmax(*$self,b);}
  
  // Binary operations with all right hand sides
  binops(const MX&)
  binops(double)
  binops(const Matrix<double> &)
  #undef binops

// all unary operations
MX __neg__(){ return - *$self;}
MX exp(){ return std::exp(*$self);}
MX log(){ return std::log(*$self);}
MX sqrt(){ return std::sqrt(*$self);}
MX sin(){ return std::sin(*$self);}
MX cos(){ return std::cos(*$self);}
MX tan(){ return std::tan(*$self);}
MX arcsin(){ return std::asin(*$self);}
MX arccos(){ return std::acos(*$self);}
MX arctan(){ return std::atan(*$self);}
MX floor(){ return std::floor(*$self);}
MX ceil(){ return std::ceil(*$self);}
MX erf(){ return std::erf(*$self);}

}

} // namespace CasADi

// Template instantiations
%template(MXVector) std::vector<CasADi::MX>;

#endif // SWIG

#endif // MX_HPP
