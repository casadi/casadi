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

namespace CasADi{

// Forward declaration
class CRSSparsityNode;
  
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

    /// Get the number of rows
    int size1() const;
    
    /// Get the number of columns
    int size2() const;

    /// Get the number of elements
    int numel() const;
    
    /// Get the number of (structural) non-zeros
    int size() const;

    /// Number of non-zeros in the upper triangular half
    int sizeU() const;

    /// Number of non-zeros in the lower triangular half
    int sizeL() const;

    /// Get a const reference to the columns of all non-zero element
    const std::vector<int>& col() const;
    
    /// Get a const reference to the rowindex of all row element
    const std::vector<int>& rowind() const;
    
    /// Get a reference to the columns of all non-zero element (copy if not unique!)
    std::vector<int>& col();
    
    /// Get a reference to the rowindex of all row element (copy if not unique!)
    std::vector<int>& rowind();
    
    /// Get the column of a non-zero element
    int col(int el) const;
    
    /// Get the index of the first non-zero element a row
    int rowind(int row) const;

    /// Resize
    void resize(int nrow, int ncol);
    
    /// Get the index of a non-zero element (copy object if necessary)
    int getNZ(int i, int j);
    
    /// Get the index of a non-zero element (return -1 if not exists)
    int getNZ(int i, int j) const;
    
    /// Get the row for each non-zero entry
    std::vector<int> getRow() const;

    /// Get the sparsity in CRS format
    void getSparsityCRS(std::vector<int>& rowind, std::vector<int> &col) const;

    /// Get the sparsity in sparse triplet format
    void getSparsity(std::vector<int>& row, std::vector<int> &col) const;
    
    /// Scalar expression
//    static const CRSSparsity scalar;

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
