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


#include "constant_mx.hpp"
#include <vector>
#include <algorithm>
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"

namespace casadi {

  ConstantMX::ConstantMX(const Sparsity& sp) {
    set_sparsity(sp);
  }

  ConstantMX::~ConstantMX() {
  }

  bool ConstantMX::is_valid_input() const {
    return nnz()==0;
  }

  casadi_int ConstantMX::n_primitives() const {
    if (nnz()==0) {
      return 0;
    } else {
      return MXNode::n_primitives();
    }
  }

  void ConstantMX::primitives(std::vector<MX>::iterator& it) const {
    if (nnz()!=0) {
      MXNode::primitives(it);
    }
  }

  template<typename T>
  void ConstantMX::split_primitives_gen(const T& x, typename std::vector<T>::iterator& it) const {
    if (nnz()!=0) {
      MXNode::split_primitives(x, it);
    }
  }

  void ConstantMX::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    split_primitives_gen<MX>(x, it);
  }

  void ConstantMX::split_primitives(const SX& x, std::vector<SX>::iterator& it) const {
    split_primitives_gen<SX>(x, it);
  }

  void ConstantMX::split_primitives(const DM& x, std::vector<DM>::iterator& it) const {
    split_primitives_gen<DM>(x, it);
  }

  template<typename T>
  T ConstantMX::join_primitives_gen(typename std::vector<T>::const_iterator& it) const {
    if (nnz()==0) {
      return T(sparsity());
    } else {
      return MXNode::join_primitives(it);
    }
  }

  MX ConstantMX::join_primitives(std::vector<MX>::const_iterator& it) const {
    if (nnz()==0) {
      return shared_from_this<MX>();
    } else {
      return MXNode::join_primitives(it);
    }
  }

  SX ConstantMX::join_primitives(std::vector<SX>::const_iterator& it) const {
    return join_primitives_gen<SX>(it);
  }

  DM ConstantMX::join_primitives(std::vector<DM>::const_iterator& it) const {
    return join_primitives_gen<DM>(it);
  }

  void ConstantMX::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = shared_from_this<MX>();
  }

 void ConstantMX::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
   MX zero_sens(size1(), size2());
   for (casadi_int d=0; d<fsens.size(); ++d) {
     fsens[d][0] = zero_sens;
   }
 }

  void ConstantMX::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                           std::vector<std::vector<MX> >& asens) const {
  }

  int ConstantMX::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    std::fill_n(res[0], nnz(), 0);
    return 0;
  }

  int ConstantMX::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    std::fill_n(res[0], nnz(), 0);
    return 0;
  }

  void ConstantDM::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res,
                            const std::vector<bool>& arg_is_ref,
                            std::vector<bool>& res_is_ref) const {
    // Print the constant
    std::string ind = g.constant(x_.nonzeros());

    if (g.elide_copy(nnz())) {
      g << g.work(res[0], nnz(), true) << " = " << ind << ";\n";
      res_is_ref[0] = true;
    } else {
      // Copy the constant to the work vector
      g << g.copy(ind, nnz(), g.work(res[0], nnz(), false)) << '\n';
    }
  }

  bool ConstantMX::__nonzero__() const {
    if (numel()!=1) casadi_error("Can only determine truth value of scalar MX.");
    if (nnz()!=1) casadi_error("Can only determine truth value of dense scalar MX.");
    return !is_zero();
  }

  ConstantMX* ConstantMX::create(const Sparsity& sp, casadi_int val) {
    if (sp.is_empty(true)) {
      return ZeroByZero::getInstance();
    } else {
      switch (val) {
      case 0: return new Constant<CompiletimeConst<0> >(sp);
      case 1: return new Constant<CompiletimeConst<1> >(sp);
      case -1: return new Constant<CompiletimeConst<(-1)> >(sp);
      default: return new Constant<RuntimeConst<casadi_int> >(sp, val);
      }
    }
  }

  ConstantMX* ConstantMX::create(const Sparsity& sp, double val) {
    if (sp.is_empty(true)) {
      return ZeroByZero::getInstance();
    } else {
      casadi_int intval = static_cast<casadi_int>(val);
      if (static_cast<double>(intval)-val==0) {
        return create(sp, intval);
      } else {
        return new Constant<RuntimeConst<double> >(sp, val);
      }
    }
  }

  ConstantMX* ConstantMX::create(const Matrix<double>& val) {
    if (val.nnz()==0) {
      return create(val.sparsity(), 0);
    } else if (val.is_scalar()) {
      return create(val.sparsity(), val.scalar());
    } else {
      // Check if all values are the same
      const std::vector<double> vdata = val.nonzeros();
      double v = vdata[0];
      for (auto&& i : vdata) {
        if (i!=v) {
          // Values not all the same
          return new ConstantDM(val);
        }
      }

      // All values identical if reached this point
      return create(val.sparsity(), v);
    }
  }

  ConstantMX* ConstantMX::create(const Sparsity& sp, const std::string& fname) {
    if (sp.nnz()==0) {
      return create(sp, 0);
    } else {
      return new ConstantFile(sp, fname);
    }
  }

  ConstantMX* ConstantMX::create(const DM& x, const std::string& name) {
    return new ConstantPool(x, name);
  }

  bool ConstantDM::is_zero() const {
    return x_.is_zero();
  }

  bool ConstantDM::is_one() const {
    return x_.is_one();
  }

  bool ConstantDM::is_minus_one() const {
    return x_.is_minus_one();
  }

  bool ConstantDM::is_eye() const {
    return x_.is_eye();
  }

  // MX ConstantMX::get_mac(const MX& y) const {
  //   if (y.is_constant()) {
  //     // Constant folding
  //     DM xv = get_DM();
  //     DM yv = y->get_DM();
  //     return mul(xv, yv);
  //   } else {
  //     return MXNode::get_mac(y);
  //   }
  // }

  MX ConstantMX::get_dot(const MX& y) const {
    if (y.is_constant()) {
      // Constant folding
      DM xv = get_DM();
      DM yv = y->get_DM();
      return dot(xv, yv);
    } else {
      return MXNode::get_dot(y);
    }
  }

  bool ConstantDM::is_equal(const MXNode* node, casadi_int depth) const {
    // Check if same node
    const ConstantDM* n = dynamic_cast<const ConstantDM*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (!std::equal(x_->begin(), x_->end(), n->x_->begin())) return false;

    return true;
  }

  std::string ZeroByZero::disp(const std::vector<std::string>& arg) const {
    return "0x0";
  }

  MX ZeroByZero::get_project(const Sparsity& sp) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const {
    casadi_assert_dev(nz.empty());
    return MX::zeros(sp);
  }

  MX ZeroByZero::get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_transpose() const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_unary(casadi_int op) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::_get_binary(casadi_int op, const MX& y, bool ScX, bool ScY) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_reshape(const Sparsity& sp) const {
    casadi_assert_dev(sp.is_empty());
    return MX::zeros(sp);
  }

  void ConstantDM::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("ConstantMX::type", 'a');
  }

  void ConstantDM::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    DM v = get_DM();
    s.pack("ConstantMX::nonzeros", v.nonzeros());
  }

  ConstantDM::ConstantDM(DeserializingStream& s) : ConstantMX(s) {
    std::vector<double> v;
    s.unpack("ConstantMX::nonzeros", v);
    x_ = DM(sparsity_, v);
  }

  void ZeroByZero::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("ConstantMX::type", 'z');
  }

  void ZeroByZero::serialize_body(SerializingStream& s) const {
   // No need to serialize body at all. All info is in header.
  }

  MXNode* ConstantMX::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("ConstantMX::type", t);
    switch (t) {
      case 'a':    return new ConstantDM(s);
      case 'f':    return new ConstantFile(s);
      case 'p':    return new ConstantPool(s);
      case 'z':    return ZeroByZero::getInstance();
      case 'D':
        return new Constant<RuntimeConst<double> >(s, RuntimeConst<double>::deserialize(s));
      case 'I':
        return new Constant<RuntimeConst<casadi_int> >(s,
                      RuntimeConst<casadi_int>::deserialize(s));
      case '0':
        return new Constant<CompiletimeConst<0> >(s, CompiletimeConst<0>::deserialize(s));
      case '1':
        return new Constant<CompiletimeConst<1> >(s, CompiletimeConst<1>::deserialize(s));
      case 'm':
        return new Constant<CompiletimeConst<(-1)> >(s, CompiletimeConst<(-1)>::deserialize(s));
      default:
        casadi_error("Error deserializing");
    }
  }

  void ConstantFile::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("ConstantFile::type", 'f');
  }

  void ConstantFile::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("ConstantFile::fname", fname_);
    s.pack("ConstantFile::x", x_);
  }

  ConstantFile::ConstantFile(DeserializingStream& s) : ConstantMX(s) {
    s.unpack("ConstantFile::fname", fname_);
    s.unpack("ConstantFile::x", x_);
  }

  ConstantFile::ConstantFile(const Sparsity& sp, const std::string& fname) :
      ConstantMX(sp), fname_(fname) {
    x_.resize(sp.nnz());
    int ret = casadi_file_slurp(fname_.c_str(), nnz(), get_ptr(x_));
    if (ret==1) casadi_error("Cannot open file '" + str(fname) + "'.");
    if (ret==2) casadi_error("Failed to read a double from '" + str(fname) + "'. "
                             "Expected " + str(sp.nnz()) + " doubles.");
  }

  std::string ConstantFile::disp(const std::vector<std::string>& arg) const {
    return "from_file('"  + fname_ + "'): " + DM(sparsity(), x_, false).get_str();
  }

  double ConstantFile::to_double() const {
    casadi_error("Not defined for ConstantFile");
  }

  Matrix<double> ConstantFile::get_DM() const {
    casadi_error("Not defined for ConstantFile");
  }

  void ConstantFile::codegen_incref(CodeGenerator& g, std::set<void*>& added) const {
    g << g.file_slurp(fname_, nnz(), g.rom_double(this)) << ";\n";
  }

  void ConstantFile::add_dependency(CodeGenerator& g) const {
    g.define_rom_double(this, nnz());
  }

  void ConstantFile::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res,
                            const std::vector<bool>& arg_is_ref,
                            std::vector<bool>& res_is_ref) const {
    if (nnz()==1) {
      g << g.workel(res[0]) << " = " << g.rom_double(this) << "[0];\n";
    } else if (g.elide_copy(nnz())) {
      g << g.work(res[0], nnz(), true) << " = " << g.rom_double(this) << ";\n";
      res_is_ref[0] = true;
    } else {
      g << g.copy(g.rom_double(this), nnz(), g.work(res[0], nnz(), false)) << '\n';
    }
  }

  ConstantPool::ConstantPool(const DM& x, const std::string& name) :
      ConstantMX(x.sparsity()), name_(name), x_(x.nonzeros()) {
  }

  std::string ConstantPool::disp(const std::vector<std::string>& arg) const {
    return "constant_pool('"  + name_ + "'): " + DM(sparsity(), x_, false).get_str();
  }

  double ConstantPool::to_double() const {
    casadi_error("Not defined for ConstantPool");
  }

  Matrix<double> ConstantPool::get_DM() const {
    casadi_error("Not defined for ConstantPool");
  }

  void ConstantPool::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res,
                            const std::vector<bool>& arg_is_ref,
                            std::vector<bool>& res_is_ref) const {
    if (nnz()==1) {
      g << g.workel(res[0]) << " = " << g.pool_double(name_) << "[0];\n";
    } else if (g.elide_copy(nnz())) {
      g << g.work(res[0], nnz(), true) << " = " << g.pool_double(name_) << ";\n";
      res_is_ref[0] = true;
    } else {
      g << g.copy(g.pool_double(name_), nnz(), g.work(res[0], nnz(), false)) << '\n';
    }
  }

  void ConstantPool::add_dependency(CodeGenerator& g) const {
    g.define_pool_double(name_, x_);
    g.add_include("string.h");
  }

  void ConstantPool::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("ConstantPool::name", name_);
    s.pack("ConstantPool::x", x_);
  }

  void ConstantPool::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("ConstantPool::type", 'p');
  }

  ConstantPool::ConstantPool(DeserializingStream& s) : ConstantMX(s) {
    s.unpack("ConstantPool::name", name_);
    s.unpack("ConstantPool::x", x_);
  }

} // namespace casadi
