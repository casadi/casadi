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

#ifndef BINARY_MX_HPP
#define BINARY_MX_HPP

#include "mx_node.hpp"

namespace CasADi{
  /** \brief Represents any binary operation that involves two matrices 
      \author Joel Andersson 
      \date 2010	
  */
  class BinaryMX : public MXNode{
  public:
    /** \brief  Constructor */
    BinaryMX(Operation op, const MX& x, const MX& y);
    
    /** \brief  Destructor */
    virtual ~BinaryMX()=0;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const{ return op_;}
    
    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Evaluate the function symbolically (SX) */
    //virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const{ return 2;}

    /** \brief Generate code for the operation (generic) */
    void generateOperationGen(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen,
			      bool el0_scalar, bool el1_scalar) const;

    //! \brief Operation
    Operation op_;
    
  };

  /// A sparse matrix-matrix binary operation
  class SparseSparseOp : public BinaryMX{
  public:
    
    /** \brief  Constructor */
    SparseSparseOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~SparseSparseOp(){};

    /** \brief  Clone function */
    virtual SparseSparseOp * clone() const{ return new SparseSparseOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    //! \brief Which argument for each nonzero
    std::vector<unsigned char> mapping_;
  };

  /// A matrix-scalar binary operation where one loops only over nonzeros of the matrix
  class NonzerosScalarOp : public BinaryMX{
  public:
    
    /** \brief  Constructor */
    NonzerosScalarOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~NonzerosScalarOp(){};

    /** \brief  Clone function */
    virtual NonzerosScalarOp * clone() const{ return new NonzerosScalarOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;
  };

  /// A scalar-matrix binary operation where one loops only over nonzeros of the matrix
  class ScalarNonzerosOp : public BinaryMX{
  public:
    
    /** \brief  Constructor */
    ScalarNonzerosOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~ScalarNonzerosOp(){};

    /** \brief  Clone function */
    virtual ScalarNonzerosOp * clone() const{ return new ScalarNonzerosOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;
  };

  /// A matrix-matrix binary operation with matching nonzeros
  class NonzerosNonzerosOp : public BinaryMX{
  public:
    
    /** \brief  Constructor */
    NonzerosNonzerosOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~NonzerosNonzerosOp(){};

    /** \brief  Clone function */
    virtual NonzerosNonzerosOp * clone() const{ return new NonzerosNonzerosOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;
  };




} // namespace CasADi


#endif // BINARY_MX_HPP
