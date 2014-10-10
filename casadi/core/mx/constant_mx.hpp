/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_CONSTANT_MX_HPP
#define CASADI_CONSTANT_MX_HPP

#include "mx_node.hpp"
#include <iomanip>
#include <iostream>

/// \cond INTERNAL

namespace casadi {

/** \brief Represents an MX that is only composed of a constant.
        \author Joel Andersson
        \date 2010-2013

        A regular user is not supposed to work with this Node class.
        This user can call MX(double) directly, or even rely on implicit typecasting.
        \sa zeros , ones
*/
  class CASADI_CORE_EXPORT ConstantMX : public MXNode {
  public:
    /// Destructor
    explicit ConstantMX(const Sparsity& sp);

    /// Destructor
    virtual ~ConstantMX() = 0;

    // Creator (all values are the same integer)
    static ConstantMX* create(const Sparsity& sp, int val);

    // Creator (all values are the same floating point value)
    static ConstantMX* create(const Sparsity& sp, double val);

    // Creator (values may be different)
    static ConstantMX* create(const Matrix<double>& val);

    /** \brief  Clone function */
    virtual ConstantMX* clone() const = 0;

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                            MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                            bool output_given);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief Get the operation */
    virtual int getOp() const { return OP_CONST;}

    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const = 0;

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const = 0;

    /// Matrix multiplication
    //    virtual MX getMultiplication(const MX& y) const;

    /// Inner product
    virtual MX getInnerProd(const MX& y) const;

    /// Return truth value of an MX
    virtual bool __nonzero__() const;
  };

  /// A constant given as a DMatrix
  class CASADI_CORE_EXPORT ConstantDMatrix : public ConstantMX {
  public:

    /** \brief  Constructor */
    explicit ConstantDMatrix(const Matrix<double>& x) : ConstantMX(x.sparsity()), x_(x) {}

    /// Destructor
    virtual ~ConstantDMatrix() {}

    /** \brief  Clone function */
    virtual ConstantDMatrix* clone() const { return new ConstantDMatrix(*this);}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const {
      x_.print(stream);
    }

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp) {
      output[0]->set(x_);
      ConstantMX::evaluateD(input, output, itmp, rtmp);
    }

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp) {
      output[0]->set(SX(x_));
      ConstantMX::evaluateSX(input, output, itmp, rtmp);
    }

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                   const std::vector<std::string>& res, CodeGenerator& gen) const;

    /** \brief  Check if a particular integer value */
    virtual bool isZero() const;
    virtual bool isOne() const;
    virtual bool isMinusOne() const;
    virtual bool isIdentity() const;

    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const {return x_.toScalar();}

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const { return x_;}

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool isEqual(const MXNode* node, int depth) const;

    /** \brief  data member */
    Matrix<double> x_;
  };

  /// A zero-by-zero matrix
  class CASADI_CORE_EXPORT ZeroByZero : public ConstantMX {
  private:
    /** \brief Private constructor (singleton design pattern) */
    explicit ZeroByZero() : ConstantMX(Sparsity::sparse(0, 0)) {
      initSingleton();
    }

  public:
    /** \brief Get a pointer to the singleton */
    static ZeroByZero* getInstance() {
      static ZeroByZero instance;
      return &instance;
    }

    /// Destructor
    virtual ~ZeroByZero() {
      destroySingleton();
    }

    /** \brief  Clone function */
    virtual ZeroByZero* clone() const { return getInstance();}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp) {}

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp) {}

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                   const std::vector<std::string>& res, CodeGenerator& gen) const {}

    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const { return 0;}

    /// Get the value (only for constant nodes)
    virtual DMatrix getMatrixValue() const { return DMatrix(); }

    /// Get densification
    virtual MX getSetSparse(const Sparsity& sp) const;

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const;

    /// Assign the nonzeros of a matrix to another matrix
    virtual MX getSetNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Transpose
    virtual MX getTranspose() const;

    /// Get a unary operation
    virtual MX getUnary(int op) const;

    /// Get a binary operation operation
    virtual MX getBinary(int op, const MX& y, bool ScX, bool ScY) const;

    /// Reshape
    virtual MX getReshape(const Sparsity& sp) const;
  };

  /** \brief Constant known at runtime */
  template<typename T>
  struct RuntimeConst {
    const T value;
    RuntimeConst() {}
    RuntimeConst(T v) : value(v) {}
  };

  /** \brief  Constant known at compiletime */
  template<int v>
  struct CompiletimeConst {
    static const int value = v;
  };

  /// A constant with all entries identical
  template<typename Value>
  class CASADI_CORE_EXPORT Constant : public ConstantMX {
  public:

    /** \brief  Constructor */
    explicit Constant(const Sparsity& sp, Value v = Value()) : ConstantMX(sp), v_(v) {}

    /// Destructor
    virtual ~Constant() {}

    /** \brief  Clone function */
    virtual Constant* clone() const { return new Constant<Value>(*this);}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp);

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                   const std::vector<std::string>& res, CodeGenerator& gen) const;

    /** \brief  Check if a particular integer value */
    virtual bool isZero() const { return v_.value==0;}
    virtual bool isOne() const { return v_.value==1;}
    virtual bool isIdentity() const { return v_.value==1 && sparsity().isDiagonal();}
    virtual bool isValue(double val) const { return v_.value==val;}

    /// Get the value (only for scalar constant nodes)
    virtual double getValue() const {
      return v_.value;
    }

    /// Get the value (only for constant nodes)
    virtual Matrix<double> getMatrixValue() const {
      return Matrix<double>(sparsity(), v_.value);
    }

    /// Get densification
    virtual MX getSetSparse(const Sparsity& sp) const;

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const;

    /// Assign the nonzeros of a matrix to another matrix
    virtual MX getSetNonzeros(const MX& y, const std::vector<int>& nz) const;

    /// Transpose
    virtual MX getTranspose() const;

    /// Get a unary operation
    virtual MX getUnary(int op) const;

    /// Get a binary operation operation
    virtual MX getBinary(int op, const MX& y, bool ScX, bool ScY) const;

    /// Reshape
    virtual MX getReshape(const Sparsity& sp) const;

    /// Create a horizontal concatenation node
    virtual MX getHorzcat(const std::vector<MX>& x) const;

    /// Create a vertical concatenation node (vectors only)
    virtual MX getVertcat(const std::vector<MX>& x) const;

    /** \brief The actual numerical value */
    Value v_;
  };

  template<typename Value>
  MX Constant<Value>::getHorzcat(const std::vector<MX>& x) const {
    // Check if all arguments have the same constant value
    for (std::vector<MX>::const_iterator i=x.begin()+1; i!=x.end(); ++i) {
      if (!(*i)->isValue(v_.value)) {
        // Not all the same value, fall back to base class
        return ConstantMX::getHorzcat(x);
      }
    }

    // Assemble the sparsity pattern
    Sparsity sp = sparsity();
    for (std::vector<MX>::const_iterator i=x.begin()+1; i!=x.end(); ++i) {
      sp.appendColumns(i->sparsity());
    }
    return MX(sp, v_.value);
  }

  template<typename Value>
  MX Constant<Value>::getVertcat(const std::vector<MX>& x) const {
    // Check if all arguments have the same constant value
    for (std::vector<MX>::const_iterator i=x.begin()+1; i!=x.end(); ++i) {
      if (!(*i)->isValue(v_.value)) {
        // Not all the same value, fall back to base class
        return ConstantMX::getVertcat(x);
      }
    }

    // Assemble the sparsity pattern
    Sparsity sp = sparsity();
    for (std::vector<MX>::const_iterator i=x.begin()+1; i!=x.end(); ++i) {
      sp.append(i->sparsity());
    }
    return MX(sp, v_.value);
  }

  template<typename Value>
  MX Constant<Value>::getReshape(const Sparsity& sp) const {
    return MX::create(new Constant<Value>(sp, v_));
  }

  template<typename Value>
  MX Constant<Value>::getTranspose() const {
    return MX::create(new Constant<Value>(sparsity().transpose(), v_));
  }

  template<typename Value>
  MX Constant<Value>::getUnary(int op) const {
    // Constant folding
    double ret(0);
    casadi_math<double>::fun(op, v_.value, 0.0, ret);
    if (operation_checker<F0XChecker>(op) || sparsity().isDense()) {
      return MX(sparsity(), ret);
    } else {
      if (v_.value==0) {
        if (isZero() && operation_checker<F0XChecker>(op)) {
          return MX(sparsity(), ret);
        } else {
          return MX::repmat(ret, size1(), size2());
        }
      }
      double ret2;
      casadi_math<double>::fun(op, 0, 0.0, ret2);
      return DMatrix(sparsity(), ret)+DMatrix(sparsity().patternInverse(), ret2);
    }
  }

  template<typename Value>
  MX Constant<Value>::getBinary(int op, const MX& y, bool ScX, bool ScY) const {
    casadi_assert(sparsity()==y.sparsity() || ScX || ScY);

    if (ScX && !operation_checker<Function0Checker>(op)) {
      double ret;
      casadi_math<double>::fun(op, size()> 0 ? v_.value: 0, 0, ret);

      if (ret!=0) {
        Sparsity f = Sparsity::dense(y.size1(), y.size2());
        MX yy = y.setSparse(f);
        return MX(f, shared_from_this<MX>())->getBinary(op, yy, false, false);
      }
    } else if (ScY && !operation_checker<F0XChecker>(op)) {
      bool grow = true;
      if (y->getOp()==OP_CONST && dynamic_cast<const ConstantDMatrix*>(y.get())==0) {
        double ret;
        casadi_math<double>::fun(op, 0, y.size()>0 ? y->getValue() : 0, ret);
        grow = ret!=0;
      }
      if (grow) {
        Sparsity f = Sparsity::dense(size1(), size2());
        MX xx = shared_from_this<MX>().setSparse(f);
        return xx->getBinary(op, MX(f, y), false, false);
      }
    }

    switch (op) {
    case OP_ADD:
      if (v_.value==0) return ScY && !y->isZero() ? MX::repmat(y, size1(), size2()) : y;
      break;
    case OP_SUB:
      if (v_.value==0) return ScY && !y->isZero() ? MX::repmat(-y, size1(), size2()) : -y;
      break;
    case OP_MUL:
      if (v_.value==1) return y;
      if (v_.value==-1) return -y;
      if (v_.value==2) return y->getUnary(OP_TWICE);
      break;
    case OP_DIV:
      if (v_.value==1) return y->getUnary(OP_INV);
      if (v_.value==-1) return -y->getUnary(OP_INV);
      break;
    case OP_POW:
      if (v_.value==0) return MX::zeros(y.sparsity());
      if (v_.value==1) return MX::ones(y.sparsity());
      if (v_.value==std::exp(1.0)) return y->getUnary(OP_EXP);
      break;
    default: break; //no rule
    }

    // Constant folding
    // NOTE: ugly, should use a function instead of a cast
    if (y->getOp()==OP_CONST && dynamic_cast<const ConstantDMatrix*>(y.get())==0) {
      double y_value = y.size()>0 ? y->getValue() : 0;
      double ret;
      casadi_math<double>::fun(op, size()> 0 ? v_.value: 0, y_value, ret);

      return MX(y.sparsity(), ret);
    }

    // Fallback
    return MXNode::getBinary(op, y, ScX, ScY);
  }

  template<typename Value>
  void Constant<Value>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                                  std::vector<int>& itmp, std::vector<double>& rtmp) {
    output[0]->set(static_cast<double>(v_.value));
    ConstantMX::evaluateD(input, output, itmp, rtmp);
  }

  template<typename Value>
  void Constant<Value>::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                                   std::vector<SXElement>& rtmp) {
    output[0]->set(SXElement(v_.value));
    ConstantMX::evaluateSX(input, output, itmp, rtmp);
  }

  template<typename Value>
  void Constant<Value>::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                          const std::vector<std::string>& res,
                                          CodeGenerator& gen) const {
    // Copy the constant to the work vector
    stream << "  for (i=0; i<" << sparsity().size() << "; ++i) ";
    stream << res.at(0) << "[i]=";
    std::ios_base::fmtflags fmtfl = stream.flags(); // get current format flags
    // full precision NOTE: hex better?
    stream << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    stream << v_.value << ";" << std::endl;
    stream.flags(fmtfl); // reset current format flags
  }

  template<typename Value>
  MX Constant<Value>::getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const {
    if (v_.value!=0) {
      // Check if any "holes"
      for (std::vector<int>::const_iterator k=nz.begin(); k!=nz.end(); ++k) {
        if (*k<0) {
          // Do not simplify
          return MXNode::getGetNonzeros(sp, nz);
        }
      }
    }
    return MX::create(new Constant<Value>(sp, v_));
  }

  template<typename Value>
  MX Constant<Value>::getSetNonzeros(const MX& y, const std::vector<int>& nz) const {
    if (y.isConstant() && y->isZero() && v_.value==0) {
      return y;
    }

    // Fall-back
    return MXNode::getSetNonzeros(y, nz);
  }

  template<typename Value>
  MX Constant<Value>::getSetSparse(const Sparsity& sp) const {
    if (isZero()) {
      return MX::create(new Constant<Value>(sp, v_));
    } else if (sp.isDense()) {
      DMatrix v = getMatrixValue();
      v.densify();
      return v;
    } else {
      return MXNode::getSetSparse(sp);
    }
  }

  template<typename Value>
  void Constant<Value>::printPart(std::ostream &stream, int part) const {
    if (sparsity().isScalar()) {
      // Print scalar
      if (sparsity().size()==0) {
        stream << "00";
      } else {
        stream << v_.value;
      }
    } else if (sparsity().isEmpty()) {
      // Print empty
      sparsity().printCompact(stream);
    } else {
      // Print value
      if (v_.value==0) {
        stream << "zeros(";
      } else if (v_.value==1) {
        stream << "ones(";
      } else if (v_.value!=v_.value) {
        stream << "nan(";
      } else if (v_.value==std::numeric_limits<double>::infinity()) {
        stream << "inf(";
      } else if (v_.value==-std::numeric_limits<double>::infinity()) {
        stream << "-inf(";
      } else {
        stream << "all_" << v_.value << "(";
      }

      // Print sparsity
      sparsity().printCompact(stream);
      stream << ")";
    }
  }

} // namespace casadi
/// \endcond


#endif // CASADI_CONSTANT_MX_HPP
