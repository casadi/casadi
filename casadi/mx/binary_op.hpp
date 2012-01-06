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

#ifndef MATRIX_MATRIX_OP_HPP
#define MATRIX_MATRIX_OP_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents any binary operation that involves two matrices 
  \author Joel Andersson 
  \date 2010	
*/
class BinaryOp : public MXNode{
  public:
    /** \brief  Constructor */
    BinaryOp(Operation op, const MX& x, const MX& y);
    
    /** \brief  Destructor */
    virtual ~BinaryOp()=0;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /// Is it a certain operation
    virtual bool isOperation(int op) const{ return op==op_;};

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Evaluate the function symbolically (SX) */
    //virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    //! \brief Operation
    Operation op_;
    
};

/// A sparse matrix-matrix binary operation
class SparseSparseOp : public BinaryOp{
  public:
    
    /** \brief  Constructor */
    SparseSparseOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~SparseSparseOp(){};

    /** \brief  Clone function */
    virtual SparseSparseOp * clone() const{ return new SparseSparseOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    //! \brief Which argument for each nonzero
    std::vector<unsigned char> mapping_;
};

/// A matrix-scalar binary operation where one loops only over nonzeros of the matrix
class NonzerosScalarOp : public BinaryOp{
  public:
    
    /** \brief  Constructor */
    NonzerosScalarOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~NonzerosScalarOp(){};

    /** \brief  Clone function */
    virtual NonzerosScalarOp * clone() const{ return new NonzerosScalarOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);
};

/// A scalar-matrix binary operation where one loops only over nonzeros of the matrix
class ScalarNonzerosOp : public BinaryOp{
  public:
    
    /** \brief  Constructor */
    ScalarNonzerosOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~ScalarNonzerosOp(){};

    /** \brief  Clone function */
    virtual ScalarNonzerosOp * clone() const{ return new ScalarNonzerosOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);
};

/// A matrix-matrix binary operation with matching nonzeros
class NonzerosNonzerosOp : public BinaryOp{
  public:
    
    /** \brief  Constructor */
    NonzerosNonzerosOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~NonzerosNonzerosOp(){};

    /** \brief  Clone function */
    virtual NonzerosNonzerosOp * clone() const{ return new NonzerosNonzerosOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);
};




} // namespace CasADi


#endif // MATRIX_MATRIX_OP_HPP
