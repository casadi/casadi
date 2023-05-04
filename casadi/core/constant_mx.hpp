/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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
#include "serializing_stream.hpp"

/// \cond INTERNAL

namespace casadi {

/** \brief Represents an MX that is only composed of a constant.

        \author Joel Andersson
        \date 2010-2013

        A regular user is not supposed to work with this Node class.
        This user can call MX(double) directly, or even rely on implicit typecasting.
        \sa zeros , ones

    \identifier{yu} */
  class CASADI_EXPORT ConstantMX : public MXNode {
  public:
    /// Destructor
    explicit ConstantMX(const Sparsity& sp);

    /// Destructor
    ~ConstantMX() override = 0;

    // Creator (all values are the same integer)
    static ConstantMX* create(const Sparsity& sp, casadi_int val);
    static ConstantMX* create(const Sparsity& sp, int val) {
      return create(sp, static_cast<casadi_int>(val));
    }

    // Creator (all values are the same floating point value)
    static ConstantMX* create(const Sparsity& sp, double val);

    // Creator (values may be different)
    static ConstantMX* create(const Matrix<double>& val);

    // Creator (values may be different)
    static ConstantMX* create(const Sparsity& sp, const std::string& fname);

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override = 0;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res,
                         casadi_int* iw, SXElem* w) const override = 0;

    /** \brief  Evaluate symbolically (MX)

        \identifier{yv} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{yw} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{yx} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward

        \identifier{yy} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{yz} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Get the operation

        \identifier{z0} */
    casadi_int op() const override { return OP_CONST;}

    /// Get the value (only for scalar constant nodes)
    double to_double() const override = 0;

    /// Get the value (only for constant nodes)
    Matrix<double> get_DM() const override = 0;

    /// Matrix multiplication
    //    virtual MX get_mac(const MX& y) const;

    /// Inner product
    MX get_dot(const MX& y) const override;

    /// Return truth value of an MX
    bool __nonzero__() const override;

    /** \brief  Check if valid function input

        \identifier{z1} */
    bool is_valid_input() const override;

    /** \brief Get the number of symbolic primitives

        \identifier{z2} */
    casadi_int n_primitives() const override;

    /** \brief Get symbolic primitives

        \identifier{z3} */
    void primitives(std::vector<MX>::iterator& it) const override;

    /** \brief Split up an expression along symbolic primitives

        \identifier{z4} */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives

        \identifier{z5} */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Detect duplicate symbolic expressions

        \identifier{z6} */
    bool has_duplicates() const override { return false;}

    /** \brief Reset the marker for an input expression

        \identifier{z7} */
    void reset_input() const override {}

    /** \brief Deserialize with type disambiguation

        \identifier{z8} */
    static MXNode* deserialize(DeserializingStream& s);

    /** \brief Deserializing constructor

        \identifier{z9} */
    explicit ConstantMX(DeserializingStream& s) : MXNode(s) {}
  };

  /// A constant given as a DM
  class CASADI_EXPORT ConstantDM : public ConstantMX {
  public:

    /** \brief  Constructor

        \identifier{za} */
    explicit ConstantDM(const Matrix<double>& x) : ConstantMX(x.sparsity()), x_(x) {}

    /// Destructor
    ~ConstantDM() override {}

    /** \brief  Print expression

        \identifier{zb} */
    std::string disp(const std::vector<std::string>& arg) const override {
      return x_.get_str();
    }

    /** \brief  Evaluate the function numerically

        \identifier{zc} */
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override {
      std::copy(x_->begin(), x_->end(), res[0]);
      return 0;
    }

    /** \brief  Evaluate the function symbolically (SX)

        \identifier{zd} */
    int eval_sx(const SXElem** arg, SXElem** res,
                         casadi_int* iw, SXElem* w) const override {
      std::copy(x_->begin(), x_->end(), res[0]);
      return 0;
    }

    /** \brief Generate code for the operation

        \identifier{ze} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief  Check if a particular integer value

        \identifier{zf} */
    bool is_zero() const override;
    bool is_one() const override;
    bool is_minus_one() const override;
    bool is_eye() const override;

    /// Get the value (only for scalar constant nodes)
    double to_double() const override {return x_.scalar();}

    /// Get the value (only for constant nodes)
    Matrix<double> get_DM() const override { return x_;}

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{zg} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** \brief  data member

        \identifier{zh} */
    Matrix<double> x_;

    /** \brief Serialize an object without type information

        \identifier{zi} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{zj} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{zk} */
    explicit ConstantDM(DeserializingStream& s);
  };

  /// A constant to be read from a file
  class CASADI_EXPORT ConstantFile : public ConstantMX {
  public:

    /** \brief  Constructor

        \identifier{zl} */
    explicit ConstantFile(const Sparsity& x, const std::string& fname);

    /// Destructor
    ~ConstantFile() override {}

    /** \brief Codegen incref

        \identifier{zm} */
    void codegen_incref(CodeGenerator& g, std::set<void*>& added) const override;

    /** \brief  Print expression

        \identifier{zn} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /// Get the value (only for scalar constant nodes)
    double to_double() const override;

    /// Get the value (only for constant nodes)
    Matrix<double> get_DM() const override;

    /** \brief  Evaluate the function numerically

        \identifier{zo} */
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override {
      std::copy(x_.begin(), x_.end(), res[0]);
      return 0;
    }

    /** \brief  Evaluate the function symbolically (SX)

        \identifier{zp} */
    int eval_sx(const SXElem** arg, SXElem** res,
                         casadi_int* iw, SXElem* w) const override {
      std::copy(x_.begin(), x_.end(), res[0]);
      return 0;
    }

    /** \brief Generate code for the operation

        \identifier{zq} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Add a dependent function

        \identifier{zr} */
    void add_dependency(CodeGenerator& g) const override;

    /** \brief file to read from

        \identifier{zs} */
    std::string fname_;

    /** \brief nonzeros

        \identifier{zt} */
    std::vector<double> x_;

    /** \brief Serialize an object without type information

        \identifier{zu} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{zv} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{zw} */
    explicit ConstantFile(DeserializingStream& s);
  };

  /// A zero-by-zero matrix
  class CASADI_EXPORT ZeroByZero : public ConstantMX {
  private:
    /** \brief Private constructor (singleton design pattern)

        \identifier{zx} */
    explicit ZeroByZero() : ConstantMX(Sparsity(0, 0)) {
      initSingleton();
    }

  public:
    /** \brief Get a pointer to the singleton

        \identifier{zy} */
    static ZeroByZero* getInstance() {
      static ZeroByZero instance;
      return &instance;
    }

    /// Destructor
    ~ZeroByZero() override {
      destroySingleton();
    }

    /** \brief  Print expression

        \identifier{zz} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate the function numerically

        \identifier{100} */
    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override {
      return 0;
    }

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res,
                         casadi_int* iw, SXElem* w) const override {
      return 0;
    }

    /** \brief Generate code for the operation

        \identifier{101} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override {}

    /// Get the value (only for scalar constant nodes)
    double to_double() const override { return 0;}

    /// Get the value (only for constant nodes)
    DM get_DM() const override { return DM(); }

    /// Get densification
    MX get_project(const Sparsity& sp) const override;

    /// Get the nonzeros of matrix
    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const override;

    /// Assign the nonzeros of a matrix to another matrix
    MX get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const override;

    /// Transpose
    MX get_transpose() const override;

    /// Get a unary operation
    MX get_unary(casadi_int op) const override;

    /// Get a binary operation operation
    MX _get_binary(casadi_int op, const MX& y, bool ScX, bool ScY) const override;

    /// Reshape
    MX get_reshape(const Sparsity& sp) const override;

    /** \brief  Check if valid function input

        \identifier{102} */
    bool is_valid_input() const override { return true;}

    /** \brief  Get the name

        \identifier{103} */
    const std::string& name() const override {
      static std::string dummyname;
      return dummyname;
    }

    /** \brief Serialize specific part of node

        \identifier{104} */
    void serialize_type(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{105} */
    void serialize_body(SerializingStream& s) const override;

  };

  /** \brief Constant known at runtime

      \identifier{106} */
  template<typename T>
  struct RuntimeConst {
    const T value;
    RuntimeConst() {}
    RuntimeConst(T v) : value(v) {}
    static char type_char();
    void serialize_type(SerializingStream& s) const {
      s.pack("Constant::value", value);
    }
    static RuntimeConst deserialize(DeserializingStream& s) {
      T v;
      s.unpack("Constant::value", v);
      return RuntimeConst(v);
    }
  };

  template<typename T>
  inline char RuntimeConst<T>::type_char() { return 'u'; }

  template<>
  inline char RuntimeConst<casadi_int>::type_char() { return 'I'; }

  template<>
  inline char RuntimeConst<double>::type_char() { return 'D'; }

  template<int v>
  struct CompiletimeConst {
    static const int value = v;
    static char type_char();
    void serialize_type(SerializingStream& s) const {}
    static CompiletimeConst deserialize(DeserializingStream& s) {
      return CompiletimeConst();
    }
  };

  template<int v>
  inline char CompiletimeConst<v>::type_char() { return 'u'; }

  template<>
  inline char CompiletimeConst<0>::type_char() { return '0'; }
  template<>
  inline char CompiletimeConst<(-1)>::type_char() { return 'm'; }
  template<>
  inline char CompiletimeConst<1>::type_char() { return '1'; }

  /// A constant with all entries identical
  template<typename Value>
  class CASADI_EXPORT Constant : public ConstantMX {
  public:

    /** \brief  Constructor

        \identifier{107} */
    explicit Constant(const Sparsity& sp, Value v = Value()) : ConstantMX(sp), v_(v) {}

    /** \brief Deserializing constructor

        \identifier{108} */
    explicit Constant(DeserializingStream& s, const Value& v);

    /// Destructor
    ~Constant() override {}

    /** \brief  Print expression

        \identifier{109} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate the function numerically

        \identifier{10a} */
    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief Generate code for the operation

        \identifier{10b} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief  Check if a particular integer value

        \identifier{10c} */
    bool is_zero() const override { return v_.value==0;}
    bool is_one() const override { return v_.value==1;}
    bool is_eye() const override { return v_.value==1 && sparsity().is_diag();}
    bool is_value(double val) const override { return v_.value==val;}

    /// Get the value (only for scalar constant nodes)
    double to_double() const override {
      return static_cast<double>(v_.value);
    }

    /// Get the value (only for constant nodes)
    Matrix<double> get_DM() const override {
      return Matrix<double>(sparsity(), to_double(), false);
    }

    /// Get densification
    MX get_project(const Sparsity& sp) const override;

    /// Get the nonzeros of matrix
    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const override;

    /// Assign the nonzeros of a matrix to another matrix
    MX get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const override;

    /// Transpose
    MX get_transpose() const override;

    /// Get a unary operation
    MX get_unary(casadi_int op) const override;

    /// Get a binary operation operation
    MX _get_binary(casadi_int op, const MX& y, bool ScX, bool ScY) const override;

    /// Reshape
    MX get_reshape(const Sparsity& sp) const override;

    /// Create a horizontal concatenation node
    MX get_horzcat(const std::vector<MX>& x) const override;

    /// Create a vertical concatenation node (vectors only)
    MX get_vertcat(const std::vector<MX>& x) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{10d} */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** \brief Serialize an object without type information

        \identifier{10e} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{10f} */
    void serialize_type(SerializingStream& s) const override;

    Value v_;
  };

  template<typename Value>
  void Constant<Value>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("ConstantMX::type", Value::type_char());
    v_.serialize_type(s);
  }

  template<typename Value>
  void Constant<Value>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
  }

  template<typename Value>
  Constant<Value>::Constant(DeserializingStream& s, const Value& v) : ConstantMX(s), v_(v) {
  }

  template<typename Value>
  MX Constant<Value>::get_horzcat(const std::vector<MX>& x) const {
    // Check if all arguments have the same constant value
    for (auto&& i : x) {
      if (!i->is_value(to_double())) {
        // Not all the same value, fall back to base class
        return ConstantMX::get_horzcat(x);
      }
    }

    // Assemble the sparsity pattern
    std::vector<Sparsity> sp;
    for (auto&& i : x) sp.push_back(i.sparsity());
    return MX(horzcat(sp), v_.value, false);
  }

  template<typename Value>
  MX Constant<Value>::get_vertcat(const std::vector<MX>& x) const {
    // Check if all arguments have the same constant value
    for (auto&& i : x) {
      if (!i->is_value(to_double())) {
        // Not all the same value, fall back to base class
        return ConstantMX::get_vertcat(x);
      }
    }

    // Assemble the sparsity pattern
    std::vector<Sparsity> sp;
    for (auto&& i : x) sp.push_back(i.sparsity());
    return MX(vertcat(sp), v_.value, false);
  }

  template<typename Value>
  MX Constant<Value>::get_reshape(const Sparsity& sp) const {
    return MX::create(new Constant<Value>(sp, v_));
  }

  template<typename Value>
  MX Constant<Value>::get_transpose() const {
    return MX::create(new Constant<Value>(sparsity().T(), v_));
  }

  template<typename Value>
  MX Constant<Value>::get_unary(casadi_int op) const {
    // Constant folding
    double ret(0);
    casadi_math<double>::fun(op, to_double(), 0.0, ret);
    if (operation_checker<F0XChecker>(op) || sparsity().is_dense()) {
      return MX(sparsity(), ret);
    } else {
      if (v_.value==0) {
        if (is_zero() && operation_checker<F0XChecker>(op)) {
          return MX(sparsity(), ret, false);
        } else {
          return repmat(MX(ret), size1(), size2());
        }
      }
      double ret2;
      casadi_math<double>::fun(op, 0, 0.0, ret2);
      return DM(sparsity(), ret, false)
        + DM(sparsity().pattern_inverse(), ret2, false);
    }
  }

  template<typename Value>
  MX Constant<Value>::_get_binary(casadi_int op, const MX& y, bool ScX, bool ScY) const {
    casadi_assert_dev(sparsity()==y.sparsity() || ScX || ScY);

    if (ScX && !operation_checker<FX0Checker>(op)) {
      double ret;
      casadi_math<double>::fun(op, nnz()> 0 ? to_double(): 0.0, 0, ret);

      if (ret!=0) {
        Sparsity f = Sparsity::dense(y.size1(), y.size2());
        MX yy = project(y, f);
        return MX(f, shared_from_this<MX>())->_get_binary(op, yy, false, false);
      }
    } else if (ScY && !operation_checker<F0XChecker>(op)) {
      bool grow = true;
      if (y->op()==OP_CONST && dynamic_cast<const ConstantDM*>(y.get())==nullptr) {
        double ret;
        casadi_math<double>::fun(op, 0, y.nnz()>0 ? y->to_double() : 0, ret);
        grow = ret!=0;
      }
      if (grow) {
        Sparsity f = Sparsity::dense(size1(), size2());
        MX xx = project(shared_from_this<MX>(), f);
        return xx->_get_binary(op, MX(f, y), false, false);
      }
    }

    switch (op) {
    case OP_ADD:
      if (v_.value==0) return ScY && !y->is_zero() ? repmat(y, size1(), size2()) : y;
      break;
    case OP_SUB:
      if (v_.value==0) return ScY && !y->is_zero() ? repmat(-y, size1(), size2()) : -y;
      break;
    case OP_MUL:
      if (v_.value==1) return y;
      if (v_.value==-1) return -y;
      if (v_.value==2) return y->get_unary(OP_TWICE);
      break;
    case OP_DIV:
      if (v_.value==1) return y->get_unary(OP_INV);
      if (v_.value==-1) return -y->get_unary(OP_INV);
      break;
    case OP_POW:
      if (v_.value==0) return MX::zeros(y.sparsity());
      if (v_.value==1) return MX::ones(y.sparsity());
      if (v_.value==std::exp(1.0)) return y->get_unary(OP_EXP);
      break;
    default: break; //no rule
    }

    // Constant folding
    // NOTE: ugly, should use a function instead of a cast
    if (y->op()==OP_CONST && dynamic_cast<const ConstantDM*>(y.get())==nullptr) {
      double y_value = y.nnz()>0 ? y->to_double() : 0;
      double ret;
      casadi_math<double>::fun(op, nnz()> 0.0 ? to_double(): 0, y_value, ret);

      return MX(y.sparsity(), ret, false);
    }

    // Fallback
    return MXNode::_get_binary(op, y, ScX, ScY);
  }

  template<typename Value>
  int Constant<Value>::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    std::fill(res[0], res[0]+nnz(), to_double());
    return 0;
  }

  template<typename Value>
  int Constant<Value>::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    std::fill(res[0], res[0]+nnz(), SXElem(v_.value));
    return 0;
  }

  template<typename Value>
  void Constant<Value>::generate(CodeGenerator& g,
                                 const std::vector<casadi_int>& arg,
                                 const std::vector<casadi_int>& res) const {
    if (nnz()==0) {
      // Quick return
    } else if (nnz()==1) {
      g << g.workel(res[0]) << " = " << g.constant(to_double()) << ";\n";
    } else {
      if (to_double()==0) {
        g << g.clear(g.work(res[0], nnz()), nnz()) << '\n';
      } else {
        g << g.fill(g.work(res[0], nnz()), nnz(), g.constant(to_double())) << '\n';
      }
    }
  }

  template<typename Value>
  MX Constant<Value>::get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const {
    if (v_.value!=0) {
      // Check if any "holes"
      for (std::vector<casadi_int>::const_iterator k=nz.begin(); k!=nz.end(); ++k) {
        if (*k<0) {
          // Do not simplify
          return MXNode::get_nzref(sp, nz);
        }
      }
    }
    return MX::create(new Constant<Value>(sp, v_));
  }

  template<typename Value>
  MX Constant<Value>::get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const {
    if (y.is_constant() && y->is_zero() && v_.value==0) {
      return y;
    }

    // Fall-back
    return MXNode::get_nzassign(y, nz);
  }

  template<typename Value>
  MX Constant<Value>::get_project(const Sparsity& sp) const {
    if (is_zero()) {
      return MX::create(new Constant<Value>(sp, v_));
    } else if (sp.is_dense()) {
      return densify(get_DM());
    } else {
      return MXNode::get_project(sp);
    }
  }

  template<typename Value>
  std::string
  Constant<Value>::disp(const std::vector<std::string>& arg) const {
    std::stringstream ss;
    if (sparsity().is_scalar()) {
      // Print scalar
      if (sparsity().nnz()==0) {
        ss << "00";
      } else {
        ss << v_.value;
      }
    } else if (sparsity().is_empty()) {
      // Print empty
      sparsity().disp(ss);
    } else {
      // Print value
      if (v_.value==0) {
        ss << "zeros(";
      } else if (v_.value==1) {
        ss << "ones(";
      } else if (v_.value!=v_.value) {
        ss << "nan(";
      } else if (v_.value==std::numeric_limits<double>::infinity()) {
        ss << "inf(";
      } else if (v_.value==-std::numeric_limits<double>::infinity()) {
        ss << "-inf(";
      } else {
        ss << "all_" << v_.value << "(";
      }

      // Print sparsity
      sparsity().disp(ss);
      ss << ")";
    }
    return ss.str();
  }

  template<typename Value>
  bool Constant<Value>::is_equal(const MXNode* node, casadi_int depth) const {
    return node->is_value(to_double()) && sparsity()==node->sparsity();
  }


} // namespace casadi
/// \endcond


#endif // CASADI_CONSTANT_MX_HPP
