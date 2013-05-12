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
#include <iomanip>
#include <iostream>

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

    // Creator (all values are the same integer)
    static ConstantMX* create(const CRSSparsity& sp, int val);

    // Creator (all values are the same floating point value)
    static ConstantMX* create(const CRSSparsity& sp, double val);

    // Creator (values may be different)
    static ConstantMX* create(const Matrix<double>& val);

    /** \brief  Clone function */
    virtual ConstantMX* clone() const = 0;

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
    
    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const = 0;

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const = 0;

    /// Matrix multiplcation
    virtual MX getMultiplication(const MX& y) const;

    /// Inner product
    virtual MX getInnerProd(const MX& y) const;

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

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    /** \brief  Check if a particular integer value */
    virtual bool isZero() const;
    virtual bool isOne() const;
    virtual bool isMinusOne() const;
    virtual bool isIdentity() const;
  
    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const{return x_.toScalar();}

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const{ return x_;}

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool isEqual(const MXNode* node, int depth) const;

    /** \brief  data member */
    Matrix<double> x_;
  };

  /** \brief Constant known at runtime */
  template<typename T>
  struct RuntimeConst{
    const T value;
    RuntimeConst(){}
    RuntimeConst(T v) : value(v){}
  };
  
  /** \brief  Constant known at compiletime */
  template<int v>
  struct CompiletimeConst{
    static const int value = v;
  };

  template<typename Value>
  class Constant : public ConstantMX{
  public:
    
    /** \brief  Constructor */
    explicit Constant(const CRSSparsity& sp, Value v = Value()) : ConstantMX(sp), v_(v){}

    /// Destructor
    virtual ~Constant(){}

    /** \brief  Clone function */
    virtual Constant* clone() const{ return new Constant<Value>(*this);}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;
    
    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;

    /** \brief  Check if a particular integer value */
    virtual bool isZero() const{ return v_.value==0;}
    virtual bool isOne() const{ return v_.value==1;}
    virtual bool isIdentity() const{ return v_.value==1 && sparsity().diagonal();}
    virtual bool isValue(double val) const{ return v_.value==val;}

    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const{
      return v_.value;
    }

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const{
      return Matrix<double>(sparsity(),v_.value);
    }

    /// Get densification
    virtual MX getSetSparse(const CRSSparsity& sp) const;

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const CRSSparsity& sp, const std::vector<int>& nz) const;

    /// Transpose
    virtual MX getTranspose() const;

    /// Get a unary operation
    virtual MX getUnary(int op) const;

    /// Get a binary operation operation
    virtual MX getBinary(int op, const MX& y, bool ScX, bool ScY) const;

    /** \brief The actual numerical value */
    Value v_;
  };

  template<typename Value>
  MX Constant<Value>::getTranspose() const{
    return MX::create(new Constant<Value>(sparsity().transpose(),v_));    
  }

  template<typename Value>
  MX Constant<Value>::getUnary(int op) const{
    // Constant folding
    double ret;
    casadi_math<double>::fun(op,v_.value,0.0,ret);
    return ret;
  }

  template<typename Value>
  MX Constant<Value>::getBinary(int op, const MX& y, bool ScX, bool ScY) const{
    switch(op){
    case OP_ADD:
      if(v_.value==0) return y;
      break;
    case OP_SUB:
      if(v_.value==0) return -y;
      break;
    case OP_MUL:
      if(v_.value==1) return y;
      if(v_.value==-1) return -y;
      if(v_.value==2) return y->getUnary(OP_TWICE);
      break;
    case OP_DIV:
      if(v_.value==1) return y->getUnary(OP_INV);
      if(v_.value==-1) return -y->getUnary(OP_INV);
      break;
    case OP_POW:
      if(v_.value==0) return MX::zeros(y.sparsity());
      if(v_.value==1) return MX::ones(y.sparsity());
      if(v_.value==std::exp(1.0)) return y->getUnary(OP_EXP);
      break;
    default: break; //no rule
    }

    // Constant folding
    if(y->getOp()==OP_CONST && dynamic_cast<const ConstantDMatrix*>(y.get())==0){ // NOTE: ugly, should use a function instead of a cast
      double y_value = y->getValue();
      double ret;
      casadi_math<double>::fun(op,v_.value,y_value,ret);
      return MX(y.sparsity(),ret);
    }

    // Fallback
    return MXNode::getBinary(op,y,ScX,ScY);
  }

  template<typename Value>
  void Constant<Value>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    output[0]->set(double(v_.value));
    ConstantMX::evaluateD(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename Value>
  void Constant<Value>::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    output[0]->set(SX(v_.value));
    ConstantMX::evaluateSX(input,output,fwdSeed,fwdSens,adjSeed,adjSens);
  }

  template<typename Value>
  void Constant<Value>::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    // Copy the constant to the work vector
    stream << "  for(i=0; i<" << sparsity().size() << "; ++i) ";
    stream << res.at(0) << "[i]=";
    std::ios_base::fmtflags fmtfl = stream.flags(); // get current format flags
    stream << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1); // full precision NOTE: hex better?
    stream << v_.value << ";" << std::endl;
    stream.flags(fmtfl); // reset current format flags
  }
  
  template<typename Value>
  MX Constant<Value>::getGetNonzeros(const CRSSparsity& sp, const std::vector<int>& nz) const{
    if(v_.value!=0){
      // Check if any "holes"
      for(std::vector<int>::const_iterator k=nz.begin(); k!=nz.end(); ++k){
        if(*k<0){
          // Do not simplify
          return MXNode::getGetNonzeros(sp,nz);
        }
      }
    }
    return MX::create(new Constant<Value>(sp,v_));
  }  

  template<typename Value>
  MX Constant<Value>::getSetSparse(const CRSSparsity& sp) const{
    if(isZero()){
      return MX::create(new Constant<Value>(sp,v_));
    } else {
      return MXNode::getSetSparse(sp);
    }
  }

  template<typename Value>
  void Constant<Value>::printPart(std::ostream &stream, int part) const{
    stream << "Const<" << v_.value << ">(";
    if(sparsity().scalar(true)){
      stream << "scalar";
    } else {
      stream << size1() << "x" << size2() << ": ";
      if(sparsity().dense()){
        stream << "dense";
      } else if(sparsity().size()==0){
        stream << "empty";          
      } else if(sparsity().diagonal()){
        stream << "diagonal";
      } else {
        stream << double(size())/sparsity().numel() << " %";
      }        
    }
    stream << ")";
  }

} // namespace CasADi


#endif // CONSTANT_MX_HPP
