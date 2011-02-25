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

#ifndef CRS_SPARSITY_HPP
#define CRS_SPARSITY_HPP

#include "../shared_object.hpp"
#include <vector>
#include <list>

namespace CasADi{

// Forward declaration
class CRSSparsityNode;
  
/** \brief General sparsity class
 * 
 * The storage format is a compressed row storage (CRS) format.\n
 * 
  In this format, the structural non-zero elements are stored in row-major order, starting from 
  the upper left corner of the matrix and ending in the lower right corner.
  
  In addition to the dimension (size1(),size2()), (i.e. the number of rows and the number of columns
  respectively), there are also two vectors of integers:
  
  1. "rowind" [length size1()+1], which contains the index to the first non-zero element on or after
     the corresponding row. All the non-zero elements of a particular i are thus the elements with 
     index el that fulfils: rowind[i] <= el < rowind[i+1].
     
  2. "col" [same length as the number of non-zero elements, size()] The columns for each of the
     structural non-zeros.
     
  Note that with this format, it is cheap to loop over all the non-zero elements of a particular row,
  constant time per elment, but expensive to jump to access a location (i,j).
  
  If the matrix is dense, i.e. length(col) == size1()*size2(), the format reduces to standard dense
  row major format, which allows access to an arbitrary element in constant time.
  
  Since the object is reference counted (it inherits from SharedObject), severl matrices are allowed
  to share the same sparsity pattern.
  
 * \see Matrix
 *
 * \author Joel Andersson 
 * \date 2010
*/
class CRSSparsity : public SharedObject{
  public:
  
    /// Default constructor, null pointer
    CRSSparsity();
    
    /// Construct a sparsity pattern (sparse/dense)
    CRSSparsity(int nrow, int ncol, bool dense=false);

    /// Construct a sparsity pattern from vectors
    CRSSparsity(int nrow, int ncol, std::vector<int> col, std::vector<int> rowind);
    
    /// Access a member function or object
    CRSSparsityNode* operator->();

    /// Const access a member function or object
    const CRSSparsityNode* operator->() const;
  
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Check if two sparsity patterns are the same
    bool operator==(const CRSSparsity& y) const;
    
    /// \name Size and element counting
    /// @{
    
    /// Get the number of rows
    int size1() const;
    
    /// Get the number of columns
    int size2() const;

    /** \brief The total number of elements, including structural zeros, i.e. size1()*size2()
        \see size()  */
    int numel() const;
    
    /** \brief Get the number of (structural) non-zeros
        \see numel() */
    int size() const;

    /** \brief Number of non-zeros in the upper triangular half, i.e. the number of elements (i,j) with j>=i */
    int sizeU() const;

    /** \brief Number of non-zeros in the lower triangular half, i.e. the number of elements (i,j) with j<=i */
    int sizeL() const;
    /// @}

    /** \brief Get a reference to col-vector, containing columns for all non-zero elements (see class description) */
    const std::vector<int>& col() const;
    
    /** \brief Get the column of a non-zero element */
    int col(int el) const;
    
    /** \brief Get a reference to the rowindex of all row element (see class description) */
    const std::vector<int>& rowind() const;

    /** \brief  Get a reference to the rowindex of row i (see class description) */
    int rowind(int i) const;

    /// Get a reference to the columns of all non-zero element (copy if not unique!)
    std::vector<int>& colRef();
    
    /// Get a reference to the rowindex of all row element (copy if not unique!)
    std::vector<int>& rowindRef();
    
    /** \brief Get the row for each non-zero entry
    Together with the col-vector, this vector gives the sparsity of the matrix in
    sparse triplet format, i.e. the row and column for each non-zero elements  */
    std::vector<int> getRow() const;
    
    /// Resize
    void resize(int nrow, int ncol);
    
    /** \brief Get the index of a non-zero element
         Add the element if it does not exist and copy object if it's not unique */
    int getNZ(int i, int j);
    
    /** \brief Get the index of an existing non-zero element
         return -1 if the element does not exists */
    int getNZ(int i, int j) const;

    /** \brief Get a set of non-zero element
         return -1 if the element does not exists */
    std::vector<int> getNZ(std::vector<int> ii, std::vector<int> jj) const;
    std::vector<int> getNZNew(std::vector<int> i, std::vector<int> j);
    std::vector<int> getNZNew(std::vector<int> i, std::vector<int> j) const;

    /// Get the sparsity in CRS format
    void getSparsityCRS(std::vector<int>& rowind, std::vector<int> &col) const;

    /// Get the sparsity in sparse triplet format
    void getSparsity(std::vector<int>& row, std::vector<int> &col) const;
    
    /// Bucket sort the elements by column
    void bucketSort(std::vector<std::list<int> >& buckets, std::vector<int>& row) const;

    /// Get a submatrix
    CRSSparsity getSub(const std::vector<int>& ii, const std::vector<int>& jj, std::vector<int>& mapping) const;
    
    /// Set a submatrix
    void setSub(const CRSSparsity& sub, const std::vector<int>& ii, const std::vector<int>& jj, std::vector<int>& mapping_nz, std::vector<int>& mapping_ind);
    
    /// Remove the structural non-zeros in a submatrix 
    CRSSparsity eraseSub(const std::vector<int>& ii, const std::vector<int>& jj, std::vector<int>& mapping) const;

    /// Transpose the matrix and get the reordering of the non-zero entries, i.e. the non-zeros of the original matrix for each non-zero of the new matrix
    CRSSparsity transpose(std::vector<int>& mapping) const;

    /** \brief Union of two sparsity patterns
    Returns the new sparsity pattern as well as a mapping with the same length as the number of non-zero elements
    The value is -1 if the non-zero comes from the first (i.e. this) object, 1 if it is from the second and 0 if from both */
    CRSSparsity patternUnion(const CRSSparsity& y, std::vector<int>& mapping) const;

    /** \brief Intersection of two sparsity patterns
    Returns the new sparsity pattern as well as a mapping with the same length as the number of non-zero elements
    The value is -1 if the non-zero comes from the first (i.e. this) object, 1 if it is from the second and 0 if from both */
    CRSSparsity patternIntersection(const CRSSparsity& y, std::vector<int>& mapping) const;

    /** \brief Sparsity pattern for a matrix-matrix product
    Returns the new sparsity pattern as well as a mapping with the same length as the number of non-zero elements
    The mapping contains a vector of the index pairs that makes up the scalar products for each non-zero */
    CRSSparsity patternProduct(const CRSSparsity& y_trans, std::vector< std::vector< std::pair<int,int> > >& mapping, bool with_mapping=true) const;

    /** \brief Sparsity pattern for a matrix-matrix product 
    No mapping */
    CRSSparsity patternProduct(const CRSSparsity& y_trans) const;
    
    /// Append another sparsity patten vertically
    void append(const CRSSparsity& sp);

    /// Reserve space
    void reserve(int nnz, int nrow);

    /// Is dense?
    bool dense() const;

    
};

#ifdef SWIG
%extend CRSSparsity {
std::string __repr__() { return $self->getRepresentation(); }
}
#endif // SWIG

#ifndef SWIG
class CRSSparsityNode : public SharedObjectNode{
  public:
    /// Construct a sparsity pattern from vectors
    CRSSparsityNode(int nrow, int ncol, std::vector<int> col, std::vector<int> rowind) : nrow_(nrow), ncol_(ncol), col_(col), rowind_(rowind){}

    /// Clone
    virtual CRSSparsityNode* clone() const{ return new CRSSparsityNode(*this); }

    /// Print representation
    virtual void repr(std::ostream &stream) const;

    /// Print description
    virtual void print(std::ostream &stream) const;

    /// Number of rows
    int nrow_;
    
    /// Number of columns
    int ncol_;
    
    /// vector of length nnz containing the columns for all the indices of the non-zero elements
    std::vector<int> col_;
    
    /// vector of length n+1 containing the index of the last non-zero element up till each row 
    std::vector<int> rowind_;
};
#endif // SWIG

} // namespace CasADi

#endif // CRS_SPARSITY_HPP
