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

    /** \brief  Create scalar constant (also implicit type conversion) */
    MX(double x);
    
    /** \brief  Create vector constant (also implicit type conversion) */
    MX(const std::vector<double> &x);
    
    /** \brief  Create sparse matrix constant (also implicit type conversion) */
    MX(const Matrix<double> &x);

    /** \brief  Matrix with all zeros */
    MX(int nrow, int ncol);
    
    /** \brief  Dense matrix filled with value val */
    MX(int nrow, int ncol, const MX& val);
    
    /** \brief  Destructor */
    ~MX();

#ifndef SWIG
    /** \brief  Get matrix non-zero */
    const MX operator[](int k) const;          
 
    /** \brief  Access matrix non-zero */
    NonZeros<MX> operator[](int k);
    
    /** \brief  Get matrix element */
    const MX operator()(int i, int j) const;
 
    /** \brief  Access matrix element */
    SubMatrix<MX> operator()(int i, int j);

    /** \brief  Get a column, or parts of a one */
    const MX operator()(const std::vector<int>& ii, int j) const;
    
    /** \brief  Access a column, or parts of a one */
    SubMatrix<MX > operator()(const std::vector<int>& ii, int j);

    /** \brief  Get a row, or parts of a one */
    const MX operator()(int i, const std::vector<int>& jj) const;
    
    /** \brief  Access a row, or parts of a one */
    SubMatrix<MX > operator()(int i, const std::vector<int>& jj);
    
    /** \brief  Get a submatrix */
    const MX operator()(const std::vector<int>& ii, const std::vector<int>& jj) const;

    /** \brief  Access a submatrix */
    SubMatrix<MX > operator()(const std::vector<int>& ii, const std::vector<int>& jj);


    
#endif // SWIG

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

    
  //@{
  /** \brief  Operators that changes the object */
    MX& operator+=(const MX &y);
    MX& operator-=(const MX &y);
    MX& operator*=(const MX &y);
    MX& operator/=(const MX &y);
  //@}

  #ifndef SWIG
  //@{
  /** \brief  Operators */
    friend MX operator+(const MX &x, const MX &y);
    friend MX operator-(const MX &x, const MX &y);
    friend MX operator*(const MX &x, const MX &y);
    friend MX operator/(const MX &x, const MX &y);
    MX operator-() const;
  //@}
  #endif // SWIG

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
  
  /** \brief  Matrix of all ones */  
  static MX ones(int nrow, int ncol);
  
  /** \brief  Identity matrix */  
  static MX eye(int n);

  
  //@{
  /// Python stuff
  const MX getitem(int k) const;
  const MX getitem(const std::vector<int>& I) const;
  const MX getitem(const std::vector< std::vector<int> > &II) const;
  MX& setitem(int k, const MX& el);
  MX& setitem(const std::vector<int> &I, const MX&  el);
  MX& setitem(const std::vector< std::vector<int> > &II, const MX&  el);
  //@}

  MX getSub(const std::vector<int>& ii, const std::vector<int>& jj) const;
  void setSub(const std::vector<int>& ii, const std::vector<int>& jj, const MX& el);
  MX getNZ(const std::vector<int>& kk) const;
  void setNZ(const std::vector<int>& kk, const MX& el);
};

//@{
/// Some typedefs
typedef double* Dptr;
typedef std::vector<double*> VDptr;
typedef std::vector<std::vector<double* > > VVDptr;
//@}


} // namespace CasADi

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

#else // SWIG

namespace CasADi {

%extend MX {
std::string __repr__() { return $self->getRepresentation();
}

// all binary operations with a particular right argument
#define binops(t) \
MX __add__(t b){    return *$self + b;} \
MX __radd__(t b){   return b + *$self;} \
MX __sub__(t b){    return *$self - b;} \
MX __rsub__(t b){   return b - *$self;} \
MX __mul__(t b){    return *$self * b;} \
MX __rmul__(t b){   return b * *$self;} \
MX __div__(t b){    return *$self / b;} \
MX __rdiv__(t b){   return b / *$self;} \
MX __pow__(t b){    return std::pow(*$self,b);} \
MX __rpow__(t b){   return std::pow(b,*$self);} \
MX fmin(t b){       return std::fmin(*$self,b);} \
MX fmax(t b){       return std::fmax(*$self,b);} \
MX prod(t b){       return prod(*$self,b);} \
MX inner_prod(t b){ return inner_prod(*$self,b);} \
MX outer_prod(t b){ return outer_prod(*$self,b);} \


// Binary operations with all right hand sides
binops(const MX&)
binops(double)
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
