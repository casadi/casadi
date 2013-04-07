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
  template<bool ScX, bool ScY>
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
    
    /** \brief Check if binary operation */
    virtual bool isBinaryOp() const { return true;}

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Evaluate the function (template)
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const{ return 2;}

    /** \brief Generate code for the operation (generic) */
    void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    //! \brief Operation
    Operation op_;
    
  };

  /// A matrix-scalar binary operation where one loops only over nonzeros of the matrix
  class MatrixScalarOp : public BinaryMX<false,true>{
  public:
    
    /** \brief  Constructor */
    MatrixScalarOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~MatrixScalarOp(){};

    /** \brief  Clone function */
    virtual MatrixScalarOp * clone() const{ return new MatrixScalarOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);
  };

  /// A scalar-matrix binary operation where one loops only over nonzeros of the matrix
  class ScalarMatrixOp : public BinaryMX<true,false>{
  public:
    
    /** \brief  Constructor */
    ScalarMatrixOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~ScalarMatrixOp(){};

    /** \brief  Clone function */
    virtual ScalarMatrixOp * clone() const{ return new ScalarMatrixOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);
  };

  /// A matrix-matrix binary operation with matching nonzeros
  class MatrixMatrixOp : public BinaryMX<false,false>{
  public:
    
    /** \brief  Constructor */
    MatrixMatrixOp(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~MatrixMatrixOp(){};

    /** \brief  Clone function */
    virtual MatrixMatrixOp * clone() const{ return new MatrixMatrixOp(*this);}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief  Evaluate the function (template) */
    template<typename T, typename MatV, typename MatVV> 
    void evaluateGen(const MatV& input, MatV& output, const MatVV& fwdSeed, MatVV& fwdSens, const MatVV& adjSeed, MatVV& adjSens);
  };




} // namespace CasADi


#endif // BINARY_MX_HPP
