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


#include "mx_node.hpp"
#include "mx_tools.hpp"
#include "../std_vector_tools.hpp"
#include <typeinfo>
#include "../matrix/matrix_tools.hpp"
#include "transpose.hpp"
#include "reshape.hpp"
#include "multiplication.hpp"
#include "subref.hpp"
#include "subassign.hpp"
#include "getnonzeros.hpp"
#include "setnonzeros.hpp"
#include "set_sparse.hpp"
#include "solve.hpp"
#include "unary_mx.hpp"
#include "binary_mx.hpp"
#include "determinant.hpp"
#include "inverse.hpp"
#include "inner_prod.hpp"
#include "norm.hpp"
#include "concat.hpp"
#include "split.hpp"
#include "assertion.hpp"

// Template implementations
#include "setnonzeros_impl.hpp"
#include "multiplication_impl.hpp"
#include "solve_impl.hpp"
#include "binary_mx_impl.hpp"

using namespace std;

namespace casadi {

  MXNode::MXNode() {
    temp = 0;
  }

  MXNode::~MXNode() {

    // Start destruction method if any of the dependencies has dependencies
    for (vector<MX>::iterator cc=dep_.begin(); cc!=dep_.end(); ++cc) {
      // Skip if constant
      if (cc->isConstant()) continue;

      // Check if there are other "owners" of the node
      if (cc->getCount()!= 1) {

        // Replace with a 0-by-0 matrix
        *cc = MX();

      } else {
        // Stack of expressions to be deleted
        std::stack<MX> deletion_stack;

        // Move the child to the deletion stack
        deletion_stack.push(*cc);
        *cc = MX();

        // Process stack
        while (!deletion_stack.empty()) {

          // Top element
          MX t = deletion_stack.top();

          // Check if the top element has dependencies with dependencies
          bool found_dep = false;

          // Start destruction method if any of the dependencies has dependencies
          for (vector<MX>::iterator ii=t->dep_.begin(); ii!=t->dep_.end(); ++ii) {

            // Skip if constant
            if (ii->isConstant()) continue;

            // Check if this is the only reference to the element
            if (ii->getCount()==1) {

              // Remove and add to stack
              deletion_stack.push(*ii);
              *ii = MX();
              found_dep = true;
              break;
            } else {
              // Replace with an element without dependencies
              *ii = MX();
            }
          }

          // Pop from stack if no dependencies found
          if (!found_dep) {
            deletion_stack.pop();
          }
        }
      }
    }
  }

  const string& MXNode::getName() const {
    throw CasadiException(string("MXNode::getName() not defined for class ")
                          + typeid(*this).name());
  }

  bool MXNode::__nonzero__() const {
    casadi_error("Can only determine truth value of a numeric MX.");

  }

  const MX& MXNode::dep(int ind) const {
    return dep_.at(ind);
  }

  MX& MXNode::dep(int ind) {
    return dep_.at(ind);
  }

  int MXNode::ndep() const {
    return dep_.size();
  }

  void MXNode::setSparsity(const Sparsity& sparsity) {
    sparsity_ = sparsity;
  }

  void MXNode::setDependencies(const MX& dep) {
    dep_.resize(1);
    dep_[0] = dep;
  }

  void MXNode::setDependencies(const MX& dep1, const MX& dep2) {
    dep_.resize(2);
    dep_[0] = dep1;
    dep_[1] = dep2;
  }

  void MXNode::setDependencies(const MX& dep1, const MX& dep2, const MX& dep3) {
    dep_.resize(3);
    dep_[0] = dep1;
    dep_[1] = dep2;
    dep_[2] = dep3;
  }

  int MXNode::addDependency(const MX& dep) {
    dep_.push_back(dep);
    return dep_.size()-1;
  }

  void MXNode::assign(const MX& d, const std::vector<int>& inz, bool add) {
    casadi_assert(0);
  }

  void MXNode::assign(const MX& d, const std::vector<int>& inz,
                      const std::vector<int>& onz, bool add) {
    casadi_assert(0);
  }

  void MXNode::setDependencies(const std::vector<MX>& dep) {
    dep_ = dep;
  }

  const Sparsity& MXNode::sparsity(int oind) const {
    casadi_assert_message(oind==0, "Index out of bounds");
    return sparsity_;
  }

  void MXNode::repr(std::ostream &stream) const {
    stream << "MX(";
    print(stream);
    stream << ")";
  }

  void MXNode::print(std::ostream &stream) const {
    long remaining_calls = MX::getMaxNumCallsInPrint();
    print(stream, remaining_calls);
  }

  void MXNode::print(std::ostream &stream, long& remaining_calls) const {
    if (remaining_calls>0) {
      remaining_calls--;
      printPart(stream, 0);
      for (int i=0; i<ndep(); ++i) {
        dep(i)->print(stream, remaining_calls);
        printPart(stream, i+1);
      }
    } else {
      stream << "...";
    }
  }

  void MXNode::printPart(std::ostream &stream, int part) const {
    casadi_assert(ndep()>1);
    casadi_assert(part>0);
    casadi_assert(part<ndep());
    stream << ", ";
  }

  Function& MXNode::getFunction() {
    throw CasadiException(string("MXNode::getFunction() not defined for class ") +
                          typeid(*this).name());
  }

  int MXNode::getFunctionOutput() const {
    throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") +
                          typeid(*this).name());
  }

  int MXNode::getFunctionInput() const {
    throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") +
                          typeid(*this).name());
  }

  void MXNode::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                         std::vector<double>& rtmp) {
    throw CasadiException(string("MXNode::evaluateD not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                          std::vector<SXElement>& rtmp) {
    throw CasadiException(string("MXNode::evaluateSX not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::evaluateMX(const MXPtrV& input, MXPtrV& output) {
    MXPtrVV fwdSeed, fwdSens, adjSeed, adjSens;
    evaluateMX(input, output, fwdSeed, fwdSens, adjSeed, adjSens, false);
  }

  void MXNode::evaluateMX(const MXPtrV& input, MXPtrV& output,
                          const MXPtrVV& fwdSeed, MXPtrVV& fwdSens,
                          const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given) {
    throw CasadiException(string("MXNode::evaluateMX not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    // By default, everything depends on everything
    bvec_t all_depend(0);

    if (fwd) {

      // Get dependencies of all inputs
      for (DMatrixPtrV::const_iterator i=input.begin(); i!=input.end(); ++i) {
        const DMatrix& m = **i;
        const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          all_depend |= v[i];
        }
      }

      // Propagate to all outputs
      for (DMatrixPtrV::iterator i=output.begin(); i!=output.end(); ++i) {
        DMatrix& m = **i;
        bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          v[i] = all_depend;
        }
      }
    } else {

      // Get dependencies of all outputs
      for (DMatrixPtrV::iterator i=output.begin(); i!=output.end(); ++i) {
        DMatrix& m = **i;
        bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          all_depend |= v[i];
          v[i] = 0;
        }
      }

      // Propagate to all inputs
      for (DMatrixPtrV::iterator i=input.begin(); i!=input.end(); ++i) {
        DMatrix& m = **i;
        bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
        for (int i=0; i<m.size(); ++i) {
          v[i] |= all_depend;
        }
      }
    }
  }

  void MXNode::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    SharedObjectNode::deepCopyMembers(already_copied);
    dep_ = deepcopy(dep_, already_copied);
  }

  MX MXNode::getOutput(int oind) const {
    casadi_assert_message(oind==0, "Output index out of bounds");
    return shared_from_this<MX>();
  }

  void MXNode::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                 const std::vector<std::string>& res, CodeGenerator& gen) const {
    stream << "#error " <<  typeid(*this).name() << ": " << arg << " => " << res << endl;
  }

  double MXNode::getValue() const {
    throw CasadiException(string("MXNode::getValue not defined for class ") + typeid(*this).name());
  }

  Matrix<double> MXNode::getMatrixValue() const {
    throw CasadiException(string("MXNode::getMatrixValue not defined for class ")
                          + typeid(*this).name());
  }

  MX MXNode::getTranspose() const {
    if (sparsity().isDense()) {
      return MX::create(new DenseTranspose(shared_from_this<MX>()));
    } else {
      return MX::create(new Transpose(shared_from_this<MX>()));
    }
  }

  MX MXNode::getReshape(const Sparsity& sp) const {
    return MX::create(new Reshape(shared_from_this<MX>(), sp));
  }


  MX MXNode::getMultiplication(const MX& y, const Sparsity& sp_z) const {
    // Get reference to transposed first argument
    MX trans_x = shared_from_this<MX>().T();

    // Form result of the right sparsity
    MX z;
    if (sp_z.isNull()) {
      Sparsity sp_z_ = sparsity().patternProductNew(y.sparsity());
      z = MX::zeros(sp_z_);
    } else {
      z = MX::zeros(sp_z);
    }
    casadi_assert_message(y.size2()==z.size2(), "Dimension error. Got y=" << y.size2()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(trans_x.size2()==z.size1(), "Dimension error. Got trans_x="
                          << trans_x.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==trans_x.size1(), "Dimension error. Got y=" << y.size1()
                          << " and trans_x" << trans_x.dimString() << ".");
    if (trans_x.isDense() && y.isDense()) {
      return MX::create(new DenseMultiplication<true, false>(z, trans_x, y));
    } else {
      return MX::create(new Multiplication<true, false>(z, trans_x, y));
    }
  }

  MX MXNode::getSolve(const MX& r, bool tr, const LinearSolver& linear_solver) const {
    if (tr) {
      return MX::create(new Solve<true>(dense(r), shared_from_this<MX>(), linear_solver));
    } else {
      return MX::create(new Solve<false>(dense(r), shared_from_this<MX>(), linear_solver));
    }
  }

  MX MXNode::getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const {
    if (nz.size()==0) {
      return MX::zeros(sp);
    } else {
      MX ret;
      if (Slice::isSlice(nz)) {
        ret = MX::create(new GetNonzerosSlice(sp, shared_from_this<MX>(), Slice(nz)));
      } else if (Slice::isSlice2(nz)) {
        Slice outer;
        Slice inner(nz, outer);
        ret = MX::create(new GetNonzerosSlice2(sp, shared_from_this<MX>(), inner, outer));
      } else {
        ret = MX::create(new GetNonzerosVector(sp, shared_from_this<MX>(), nz));
      }
      simplify(ret);
      return ret;
    }
  }

  MX MXNode::getSetNonzeros(const MX& y, const std::vector<int>& nz) const {
    // Check if any element needs to be set at all
    bool set_any = false;
    for (vector<int>::const_iterator i=nz.begin(); i!=nz.end() && !set_any; ++i) {
      set_any = *i >= 0;
    }

    // Quick return
    if (!set_any) {
      return y;
    }

    // Check if slice
    MX ret;
    if (Slice::isSlice(nz)) {
      ret = MX::create(new SetNonzerosSlice<false>(y, shared_from_this<MX>(), Slice(nz)));
    } else if (Slice::isSlice2(nz)) {
      Slice outer;
      Slice inner(nz, outer);
      ret = MX::create(new SetNonzerosSlice2<false>(y, shared_from_this<MX>(), inner, outer));
    } else {
      ret = MX::create(new SetNonzerosVector<false>(y, shared_from_this<MX>(), nz));
    }
    simplify(ret);
    return ret;
  }


  MX MXNode::getAddNonzeros(const MX& y, const std::vector<int>& nz) const {
    if (nz.size()==0 || isZero()) {
      return y;
    } else {
      MX ret;
      if (Slice::isSlice(nz)) {
        ret = MX::create(new SetNonzerosSlice<true>(y, shared_from_this<MX>(), Slice(nz)));
      } else if (Slice::isSlice2(nz)) {
        Slice outer;
        Slice inner(nz, outer);
        ret = MX::create(new SetNonzerosSlice2<true>(y, shared_from_this<MX>(), inner, outer));
      } else {
        ret = MX::create(new SetNonzerosVector<true>(y, shared_from_this<MX>(), nz));
      }
      simplify(ret);
      return ret;
    }
  }

  MX MXNode::getSetSparse(const Sparsity& sp) const {
    return MX::create(new SetSparse(shared_from_this<MX>(), sp));
  }

  MX MXNode::getSubRef(const Slice& i, const Slice& j) const {
    return MX::create(new SubRef(shared_from_this<MX>(), i, j));
  }

  MX MXNode::getSubAssign(const MX& y, const Slice& i, const Slice& j) const {
    return MX::create(new SubAssign(shared_from_this<MX>(), y, i, j));
  }

  MX MXNode::getUnary(int op) const {
    if (operation_checker<F0XChecker>(op) && isZero()) {
      // If identically zero
      return MX::zeros(sparsity());
    } else {
      // Create a new node
      return MX::create(new UnaryMX(Operation(op), shared_from_this<MX>()));
    }
  }

  MX MXNode::getBinarySwitch(int op, const MX& y) const {
    // Make sure that dimensions match
    casadi_assert_message((sparsity().isScalar(false) || y.isScalar() ||
                           (sparsity().size2()==y.size2() && size1()==y.size1())),
                          "Dimension mismatch." << "lhs is " << sparsity().dimString()
                          << ", while rhs is " << y.dimString());

    // Create binary node
    if (sparsity().isScalar(false)) {
      if (size()==0) {
        return toMatrix(MX(0)->getBinary(op, y, true, false), y.sparsity());
      } else {
        return toMatrix(getBinary(op, y, true, false), y.sparsity());
      }
    } else if (y.isScalar()) {
      if (y.size()==0) {
        return toMatrix(getBinary(op, MX(0), false, true), sparsity());
      } else {
        return toMatrix(getBinary(op, y, false, true), sparsity());
      }
    } else {
      casadi_assert_message(sparsity().shape() == y.sparsity().shape(), "Dimension mismatch.");
      if (sparsity()==y.sparsity()) {
        // Matching sparsities
        return getBinary(op, y, false, false);
      } else {
        // Get the sparsity pattern of the result
        // (ignoring structural zeros giving rise to nonzero result)
        const Sparsity& x_sp = sparsity();
        const Sparsity& y_sp = y.sparsity();
        Sparsity r_sp = x_sp.patternCombine(y_sp, operation_checker<F0XChecker>(op),
                                            operation_checker<Function0Checker>(op));

        // Project the arguments to this sparsity
        MX xx = shared_from_this<MX>().setSparse(r_sp);
        MX yy = y.setSparse(r_sp);
        return xx->getBinary(op, yy, false, false);
      }
    }
  }

  MX MXNode::getBinary(int op, const MX& y, bool scX, bool scY) const {
    casadi_assert(sparsity()==y.sparsity() || scX || scY);

    if (CasadiOptions::simplification_on_the_fly) {

      // If identically zero due to one argumebt being zero
      if ((operation_checker<F0XChecker>(op) && isZero()) ||
         (operation_checker<Function0Checker>(op) && y->isZero())) {
        return MX::zeros(sparsity());
      }

      // Handle special operations (independent of type)
      switch (op) {
      case OP_ADD:
        if (y.isEqual(this, maxDepth())) return getUnary(OP_TWICE);
        break;
      case OP_SUB:
      case OP_NE:
      case OP_LT:
        if (y.isEqual(this, maxDepth())) return MX::zeros(sparsity());
        break;
      case OP_DIV:
        if (y->isZero()) return MX::nan(sparsity());
        // fall-through
      case OP_EQ:
      case OP_LE:
        if (y.isEqual(this, maxDepth())) return MX::ones(sparsity());
        break;
      case OP_MUL:
        if (y.isEqual(this, maxDepth())) return getUnary(OP_SQ);
        break;
      default: break; // no rule
      }

      // Handle special cases for the second argument
      switch (y->getOp()) {
      case OP_CONST:
        // Make the constant the first argument, if possible
        if (getOp()!=OP_CONST && operation_checker<CommChecker>(op)) {
              return y->getBinary(op, shared_from_this<MX>(), scY, scX);
        } else {
          switch (op) {
          case OP_CONSTPOW:
            if (y->isValue(-1)) return getUnary(OP_INV);
            else if (y->isValue(0)) return MX::ones(sparsity());
            else if (y->isValue(1)) return shared_from_this<MX>();
            else if (y->isValue(2)) return getUnary(OP_SQ);
            break;
          case OP_ADD:
          case OP_SUB:
            if (y->isZero())
                return scX ? MX::repmat(shared_from_this<MX>(), y.shape()) : shared_from_this<MX>();
            break;
          case OP_MUL:
            if (y->isValue(1)) return shared_from_this<MX>();
            break;
          case OP_DIV:
            if (y->isValue(1)) return shared_from_this<MX>();
            else if (y->isValue(0.5)) return getUnary(OP_TWICE);
            break;
          default: break; // no rule
          }
        }
        break;
      case OP_NEG:
        if (op==OP_ADD) {
          return getBinary(OP_SUB, y->dep(), scX, scY);
        } else if (op==OP_SUB) {
          return getBinary(OP_ADD, y->dep(), scX, scY);
        } else if (op==OP_MUL) {
          return -getBinary(OP_MUL, y->dep(), scX, scY);
        } else if (op==OP_DIV) {
          return -getBinary(OP_DIV, y->dep(), scX, scY);
        }
        break;
      case OP_INV:
        if (op==OP_MUL) {
          return getBinary(OP_DIV, y->dep(), scX, scY);
        } else if (op==OP_DIV) {
          return getBinary(OP_MUL, y->dep(), scX, scY);
        }
        break;
      default: break; // no rule
      }

    }

    if (scX) {
      // Check if it is ok to loop over nonzeros only
      if (y.isDense() || operation_checker<Function0Checker>(op)) {
        // Loop over nonzeros
        return MX::create(new BinaryMX<true, false>(Operation(op), shared_from_this<MX>(), y));
      } else {
        // Put a densification node in between
        return getBinary(op, dense(y), true, false);
      }
    } else if (scY) {
      // Check if it is ok to loop over nonzeros only
      if (sparsity().isDense() || operation_checker<F0XChecker>(op)) {
        // Loop over nonzeros
        return MX::create(new BinaryMX<false, true>(Operation(op), shared_from_this<MX>(), y));
      } else {
        // Put a densification node in between
        return dense(shared_from_this<MX>())->getBinary(op, y, false, true);
      }
    } else {
      // Loop over nonzeros only
      MX rr = MX::create(new BinaryMX<false, false>(Operation(op), shared_from_this<MX>(), y));

      // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
      if (!rr.isDense() && !operation_checker<F00Checker>(op)) {
        // Get the value for the structural zeros
        double fcn_0(0);
        casadi_math<double>::fun(op, 0, 0, fcn_0);
        rr.densify(fcn_0);
      }
      return rr;
    }
  }

  Matrix<int> MXNode::mapping() const {
    throw CasadiException(string("MXNode::mapping not defined for class ") + typeid(*this).name());
  }

  bool MXNode::sameOpAndDeps(const MXNode* node, int depth) const {
    if (getOp()!=node->getOp() || ndep()!=node->ndep())
      return false;
    for (int i=0; i<ndep(); ++i) {
      if (!dep(i).isEqual(node->dep(i), depth-1))
        return false;
    }
    return true;
  }

  MX MXNode::getAssertion(const MX& y, const std::string &fail_message) const {
    return MX::create(new Assertion(shared_from_this<MX>(), y, fail_message));
  }

  MX MXNode::getDeterminant() const {
    return MX::create(new Determinant(shared_from_this<MX>()));
  }

  MX MXNode::getInverse() const {
    return MX::create(new Inverse(shared_from_this<MX>()));
  }


  MX MXNode::getInnerProd(const MX& y) const {
    casadi_assert_message(
      size2()==y.size2() && size1()==y.size1(),
      "MXNode::inner_prod: Dimension mismatch. inner_prod requires its "
      "two arguments to have equal shapes, but got ("
      << size2() << ", " << size1() << ") and ("
      << y.size2() << ", " << y.size1() << ").");
    if (sparsity()==y.sparsity()) {
      if (sparsity().size()==0) {
        return 0;
      } else if (sparsity().isScalar()) {
        return getBinarySwitch(OP_MUL, y);
      } else {
        return MX::create(new InnerProd(shared_from_this<MX>(), y));
      }
    } else {
      // Project to pattern intersection
      Sparsity sp = sparsity().patternIntersection(y.sparsity());
      MX xx = shared_from_this<MX>().setSparse(sp);
      MX yy = y.setSparse(sp);
      return xx->getInnerProd(yy);
    }
  }

  MX MXNode::getNormF() const {
    return MX::create(new NormF(shared_from_this<MX>()));
  }

  MX MXNode::getNorm2() const {
    return MX::create(new Norm2(shared_from_this<MX>()));
  }

  MX MXNode::getNormInf() const {
    return MX::create(new NormInf(shared_from_this<MX>()));
  }

  MX MXNode::getNorm1() const {
    return MX::create(new Norm1(shared_from_this<MX>()));
  }

  MX MXNode::getHorzcat(const std::vector<MX>& x) const {
    // Check if there is any existing horzcat operation
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (i->getOp()==OP_HORZCAT) {
        // Split up
        vector<MX> x_split(x.begin(), i);
        for (; i!=x.end(); ++i) {
          if (i->getOp()==OP_HORZCAT) {
            x_split.insert(x_split.end(), (*i)->dep_.begin(), (*i)->dep_.end());
          } else {
            x_split.push_back(*i);
          }
        }
        return horzcat(x_split);
      }
    }

    // Create a Horzcat node
    return MX::create(new Horzcat(x));
  }

  MX MXNode::getDiagcat(const std::vector<MX>& x) const {
    // Create a Horzcat node
    return MX::create(new Diagcat(x));
  }

  MX MXNode::getVertcat(const std::vector<MX>& x) const {
    // Check if there is any existing vertcat operation
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (i->getOp()==OP_VERTCAT) {
        // Split up
        vector<MX> x_split(x.begin(), i);
        for (; i!=x.end(); ++i) {
          if (i->getOp()==OP_VERTCAT) {
            x_split.insert(x_split.end(), (*i)->dep_.begin(), (*i)->dep_.end());
          } else {
            x_split.push_back(*i);
          }
        }
        return vertcat(x_split);
      }
    }

    return MX::create(new Vertcat(x));
  }

  std::vector<MX> MXNode::getHorzsplit(const std::vector<int>& output_offset) const {
    if (isZero()) {
      std::vector<MX> ret =
          MX::createMultipleOutput(new Horzsplit(shared_from_this<MX>(), output_offset));
      for (int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    std::vector<MX> ret =
        MX::createMultipleOutput(new Horzsplit(shared_from_this<MX>(), output_offset));

    if (CasadiOptions::simplification_on_the_fly) {
      // Simplify horzsplit(horzcat)
      if (getOp()==OP_HORZCAT) {
        int offset_deps = 0;
        int j = 0;
        for (int i=0;i<output_offset.size();++i) {
          while (offset_deps<output_offset[i]) { offset_deps+=dep(j).size2();++j; }
          if (j>=ndep()) j = ndep()-1;
          if (output_offset[i]==offset_deps &&
              (i+1<output_offset.size()?output_offset[i+1]:size2()) ==
               offset_deps +dep(j).size2()) {
            // Aligned with vertcat dependency
            ret[i] = dep(j);
          }
        }
      }
    }
    return ret;
  }

  std::vector<MX> MXNode::getDiagsplit(const std::vector<int>& offset1,
                                       const std::vector<int>& offset2) const {
    if (isZero()) {
      std::vector<MX> ret =
          MX::createMultipleOutput(new Diagsplit(shared_from_this<MX>(), offset1, offset2));
      for (int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    std::vector<MX> ret =
        MX::createMultipleOutput(new Diagsplit(shared_from_this<MX>(), offset1, offset2));

    return ret;
  }

  std::vector<MX> MXNode::getVertsplit(const std::vector<int>& output_offset) const {
    if (isZero()) {
      std::vector<MX> ret =
          MX::createMultipleOutput(new Vertsplit(shared_from_this<MX>(), output_offset));
      for (int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    std::vector<MX> ret =
        MX::createMultipleOutput(new Vertsplit(shared_from_this<MX>(), output_offset));

    if (CasadiOptions::simplification_on_the_fly) {
      // Simplify vertsplit(vertcat)
      if (getOp()==OP_VERTCAT) {
        int offset_deps = 0;
        int j = 0;
        for (int i=0;i<output_offset.size();++i) {
          while (offset_deps<output_offset[i]) { offset_deps+=dep(j).size1();++j; }
          if (j>=ndep()) j = ndep()-1;
          if (output_offset[i]==offset_deps &&
              (i+1<output_offset.size()?output_offset[i+1]:size1()) ==
               offset_deps +dep(j).size1()) {
            // Aligned with vertcat dependency
            ret[i] = dep(j);
          }
        }
      }
    }
    return ret;
  }

  void MXNode::clearVector(const std::vector<std::vector<MX*> > v) {
    for (int i=0; i<v.size(); ++i) {
      for (int j=0; j<v[i].size(); ++j) {
        if (v[i][j]!= 0) {
          *v[i][j] = MX();
        }
      }
    }
  }


} // namespace casadi
