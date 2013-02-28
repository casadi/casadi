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

#ifndef CONSTANT_MX_HPP
#define CONSTANT_MX_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Represents an MX that is only composed of a constant.
	\author Joel Andersson 
	\date 2010-2013

	A regular user is not supposed to work with this Node class.
	This user can call MX(double) directly, or even rely on implicit typecasting.
	\sa zeros , ones
*/
  class ConstantMX : public MXNode{
  public:
    /// Destructor
    explicit ConstantMX(const CRSSparsity& sp);

    /// Destructor
    virtual ~ConstantMX() = 0;

    /** \brief  Clone function */
    virtual ConstantMX* clone() const = 0;

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);
   
    /** \brief Get the operation */
    virtual int getOp() const{ return OP_CONST;}
    
    /// Return truth value of an MX
    virtual bool __nonzero__() const;
  };

  class ConstantDMatrix : public ConstantMX{
  public:

    /** \brief  Constructor */
    explicit ConstantDMatrix(const Matrix<double>& x) : ConstantMX(x.sparsity()), x_(x){}
    
    /// Destructor
    virtual ~ConstantDMatrix(){}

    /** \brief  Clone function */
    virtual ConstantDMatrix* clone() const{ return new ConstantDMatrix(*this);}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const{
      x_.print(stream);
    }
    
    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){ 
      output[0]->set(x_);
      ConstantMX::evaluateD(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
    }

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
      output[0]->set(SXMatrix(x_));
      ConstantMX::evaluateSX(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
    }

    /** \brief  Check if a particular integer value */
    virtual bool isZero() const;
    virtual bool isOne() const;
    virtual bool isMinusOne() const;
    virtual bool isIdentity() const;
  
    /** \brief  data member */
    Matrix<double> x_;
  };

  template<int value>
  class ConstantInt : public ConstantMX{
  public:
    
    /** \brief  Constructor */
    ConstantInt(const CRSSparsity& sp) : ConstantMX(sp){}

    /// Destructor
    virtual ~ConstantInt(){}

    /** \brief  Clone function */
    virtual ConstantInt* clone() const{ return new ConstantInt<value>(*this);}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const{
      stream << value;
    }
    
    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
      output[0]->set(double(value));
      ConstantMX::evaluateD(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
    }

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
      output[0]->set(SX(value));
      ConstantMX::evaluateSX(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
    }

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
      // Copy the constant to the work vector
      stream << "  for(i=0; i<" << sparsity().size() << "; ++i) ";
      stream << res.at(0) << "[i]=";
      stream << value << ";" << std::endl;
    }

    /** \brief  Check if a particular integer value */
    virtual bool isZero() const{ return value==0;}
    virtual bool isOne() const{ return value==1;}
    virtual bool isMinusOne() const{ return value==-1;}
    virtual bool isIdentity() const{ return value==1 && sparsity().diagonal();}
  };


} // namespace CasADi


#endif // CONSTANT_MX_HPP
