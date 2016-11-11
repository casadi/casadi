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
#include "../std_vector_tools.hpp"
#include <typeinfo>
#include "transpose.hpp"
#include "reshape.hpp"
#include "multiplication.hpp"
#include "bilin.hpp"
#include "rank1.hpp"
#include "subref.hpp"
#include "subassign.hpp"
#include "getnonzeros.hpp"
#include "setnonzeros.hpp"
#include "project.hpp"
#include "solve.hpp"
#include "unary_mx.hpp"
#include "binary_mx.hpp"
#include "determinant.hpp"
#include "inverse.hpp"
#include "dot.hpp"
#include "norm.hpp"
#include "concat.hpp"
#include "split.hpp"
#include "assertion.hpp"
#include "monitor.hpp"
#include "repmat.hpp"
#include "casadi_find.hpp"

// Template implementations
#include "setnonzeros_impl.hpp"
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
      if (cc->is_constant()) continue;

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
            if (ii->is_constant()) continue;

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

  int MXNode::n_primitives() const {
    throw CasadiException(string("MXNode::n_primitives() not defined for class ")
                          + typeid(*this).name());
  }

  bool MXNode::has_duplicates() {
    throw CasadiException(string("MXNode::has_duplicates() not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::resetInput() {
    throw CasadiException(string("MXNode::resetInput() not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::primitives(vector<MX>::iterator& it) const {
    throw CasadiException(string("MXNode::primitives() not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::split_primitives(const MX& x, vector<MX>::iterator& it) const {
    throw CasadiException(string("MXNode::split_primitives() not defined for class ")
                          + typeid(*this).name());
  }

  MX MXNode::join_primitives(vector<MX>::const_iterator& it) const {
    throw CasadiException(string("MXNode::join_primitives() not defined for class ")
                          + typeid(*this).name());
  }

  const string& MXNode::name() const {
    throw CasadiException(string("MXNode::name() not defined for class ")
                          + typeid(*this).name());
  }

  std::string MXNode::type_name() const {
    return typeid(*this).name();
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

  void MXNode::assign(const MX& d, const vector<int>& inz, bool add) {
    casadi_assert(0);
  }

  void MXNode::assign(const MX& d, const vector<int>& inz,
                      const vector<int>& onz, bool add) {
    casadi_assert(0);
  }

  void MXNode::setDependencies(const vector<MX>& dep) {
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
    // Find out which noded can be inlined
    std::map<const MXNode*, int> nodeind;
    can_inline(nodeind);

    // Print expression
    vector<string> intermed;
    string s = print_compact(nodeind, intermed);

    // Print intermediate expressions
    for (int i=0; i<intermed.size(); ++i)
      stream << "@" << (i+1) << "=" << intermed[i] << ", ";

    // Print this
    stream << s;
  }

  void MXNode::can_inline(std::map<const MXNode*, int>& nodeind) const {
    // Add or mark node in map
    std::map<const MXNode*, int>::iterator it=nodeind.find(this);
    if (it==nodeind.end()) {
      // First time encountered, mark inlined
      nodeind.insert(it, make_pair(this, 0));

      // Handle dependencies with recursion
      for (int i=0; i<ndep(); ++i) {
        dep(i)->can_inline(nodeind);
      }
    } else if (it->second==0 && op()!=OP_PARAMETER) {
      // Node encountered before, do not inline (except if symbolic primitive)
      it->second = -1;
    }
  }

  std::string MXNode::print_compact(std::map<const MXNode*, int>& nodeind,
                                   vector<std::string>& intermed) const {
    // Get reference to node index
    int& ind = nodeind[this];

    // If positive, already in intermediate expressions
    if (ind>0) return "@" + CodeGenerator::to_string(ind);

    // Get expressions for dependencies
    vector<string> arg(ndep());
    for (int i=0; i<arg.size(); ++i) {
      arg[i] = dep(i)->print_compact(nodeind, intermed);
    }

    // Get expression for this
    string s = print(arg);

    // Decide what to do with the expression
    if (ind==0) {
      // Inline expression
      return s;
    } else {
      // Add to list of intermediate expressions and return reference
      intermed.push_back(s);
      ind = intermed.size(); // For subsequent references
      return "@" + CodeGenerator::to_string(ind);
    }
  }

  const Function& MXNode::getFunction(int i) const {
    throw CasadiException(string("MXNode::getFunction() not defined for class ") +
                          typeid(*this).name());
  }

  int MXNode::getFunctionOutput() const {
    throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") +
                          typeid(*this).name());
  }

  int MXNode::getFunction_input() const {
    throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") +
                          typeid(*this).name());
  }

  void MXNode::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    throw CasadiException(string("MXNode::eval not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    throw CasadiException(string("MXNode::eval_sx not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::eval_mx(const vector<MX>& arg, vector<MX>& res) {
    throw CasadiException(string("MXNode::eval_mx not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::evalFwd(const vector<vector<MX> >& fseed,
                       vector<vector<MX> >& fsens) {
    throw CasadiException(string("MXNode::evalFwd not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::evalAdj(const vector<vector<MX> >& aseed,
                       vector<vector<MX> >& asens) {
    throw CasadiException(string("MXNode::evalAdj not defined for class ")
                          + typeid(*this).name());
  }

  void MXNode::sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    // By default, everything depends on everything
    bvec_t all_depend(0);

    // Get dependencies of all inputs
    for (int k=0; k<ndep(); ++k) {
      const bvec_t* v = arg[k];
      for (int i=0; i<dep(k).nnz(); ++i) {
        all_depend |= v[i];
      }
    }

    // Propagate to all outputs
    for (int k=0; k<nout(); ++k) {
      bvec_t* v = res[k];
      for (int i=0; i<sparsity(k).nnz(); ++i) {
        v[i] = all_depend;
      }
    }
  }

  void MXNode::sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    // By default, everything depends on everything
    bvec_t all_depend(0);

    // Get dependencies of all outputs
    for (int k=0; k<nout(); ++k) {
      bvec_t* v = res[k];
      for (int i=0; i<sparsity(k).nnz(); ++i) {
        all_depend |= v[i];
        v[i] = 0;
      }
    }

    // Propagate to all inputs
    for (int k=0; k<ndep(); ++k) {
      bvec_t* v = arg[k];
      for (int i=0; i<dep(k).nnz(); ++i) {
        v[i] |= all_depend;
      }
    }
  }

  MX MXNode::getOutput(int oind) const {
    casadi_assert_message(oind==0, "Output index out of bounds");
    return shared_from_this<MX>();
  }

  void MXNode::generate(CodeGenerator& g, const std::string& mem,
                        const vector<int>& arg, const vector<int>& res) const {
    casadi_warning("Cannot code generate MX nodes of type " + type_name() +
                   "The generation will proceed, but compilation of the code will "
                   "not be possible.");
    g.body << "#error " <<  type_name() << ": " << arg << " => " << res << endl;
  }

  double MXNode::to_double() const {
    casadi_error(string("MXNode::to_double not defined for class ") + typeid(*this).name());
  }

  Matrix<double> MXNode::getMatrixValue() const {
    casadi_error(string("MXNode::getMatrixValue not defined for class ")
                 + typeid(*this).name());
  }

  MX MXNode::getTranspose() const {
    if (sparsity().is_scalar()) {
      return shared_from_this<MX>();
    } else if (sparsity().is_vector()) {
      return getReshape(sparsity().T());
    } else if (sparsity().is_dense()) {
      return MX::create(new DenseTranspose(shared_from_this<MX>()));
    } else {
      return MX::create(new Transpose(shared_from_this<MX>()));
    }
  }

  MX MXNode::getReshape(const Sparsity& sp) const {
    casadi_assert(sp.isReshape(sparsity()));
    if (sp==sparsity()) {
      return shared_from_this<MX>();
    } else {
      return MX::create(new Reshape(shared_from_this<MX>(), sp));
    }
  }


  MX MXNode::getMultiplication(const MX& y, const MX& z) const {
    // Get reference to transposed first argument
    MX x = shared_from_this<MX>();

    casadi_assert_message(y.size2()==z.size2(), "Dimension error. Got y=" << y.size2()
                          << " and z=" << z.dim() << ".");
    casadi_assert_message(x.size1()==z.size1(), "Dimension error. Got x="
                          << x.dim() << " and z=" << z.dim() << ".");
    casadi_assert_message(y.size1()==x.size2(), "Dimension error. Got y=" << y.size1()
                          << " and x" << x.dim() << ".");
    if (x.is_dense() && y.is_dense() && z.is_dense()) {
      return MX::create(new DenseMultiplication(z, x, y));
    } else {
      return MX::create(new Multiplication(z, x, y));
    }
  }

  MX MXNode::getBilin(const MX& x, const MX& y) const {
    return MX::create(new Bilin(shared_from_this<MX>(), x, y));
  }

  MX MXNode::getRank1(const MX& alpha, const MX& x, const MX& y) const {
    return MX::create(new Rank1(shared_from_this<MX>(), alpha, x, y));
  }

  MX MXNode::getSolve(const MX& r, bool tr, const Linsol& linear_solver) const {
    if (tr) {
      return MX::create(new Solve<true>(densify(r), shared_from_this<MX>(), linear_solver));
    } else {
      return MX::create(new Solve<false>(densify(r), shared_from_this<MX>(), linear_solver));
    }
  }

  MX MXNode::getGetNonzeros(const Sparsity& sp, const vector<int>& nz) const {
    if (nz.size()==0) {
      return MX::zeros(sp);
    } else {
      MX ret;
      if (is_slice(nz)) {
        ret = MX::create(new GetNonzerosSlice(sp, shared_from_this<MX>(), to_slice(nz)));
      } else if (is_slice2(nz)) {
        pair<Slice, Slice> sl = to_slice2(nz);
        ret = MX::create(new GetNonzerosSlice2(sp, shared_from_this<MX>(), sl.first, sl.second));
      } else {
        ret = MX::create(new GetNonzerosVector(sp, shared_from_this<MX>(), nz));
      }
      return simplify(ret);
    }
  }

  MX MXNode::getSetNonzeros(const MX& y, const vector<int>& nz) const {
    // Check if any element needs to be set at all
    bool set_any = false;
    for (auto i=nz.begin(); i!=nz.end() && !set_any; ++i) {
      set_any = *i >= 0;
    }

    // Quick return
    if (!set_any) {
      return y;
    }

    // Check if slice
    MX ret;
    if (is_slice(nz)) {
      ret = MX::create(new SetNonzerosSlice<false>(y, shared_from_this<MX>(), to_slice(nz)));
    } else if (is_slice2(nz)) {
      pair<Slice, Slice> sl = to_slice2(nz);
      ret = MX::create(new SetNonzerosSlice2<false>(y, shared_from_this<MX>(),
                                                    sl.first, sl.second));
    } else {
      ret = MX::create(new SetNonzerosVector<false>(y, shared_from_this<MX>(), nz));
    }
    return simplify(ret);
  }


  MX MXNode::getAddNonzeros(const MX& y, const vector<int>& nz) const {
    if (nz.size()==0 || is_zero()) {
      return y;
    } else {
      MX ret;
      if (is_slice(nz)) {
        ret = MX::create(new SetNonzerosSlice<true>(y, shared_from_this<MX>(), to_slice(nz)));
      } else if (is_slice2(nz)) {
        pair<Slice, Slice> sl = to_slice2(nz);
        ret = MX::create(new SetNonzerosSlice2<true>(y, shared_from_this<MX>(),
                                                     sl.first, sl.second));
      } else {
        ret = MX::create(new SetNonzerosVector<true>(y, shared_from_this<MX>(), nz));
      }
      return simplify(ret);
    }
  }

  MX MXNode::getProject(const Sparsity& sp) const {
    if (sp==sparsity()) {
      return shared_from_this<MX>();
    } else if (sp.nnz()==0) {
      return MX::zeros(sp);
    } else {
      return MX::create(new Project(shared_from_this<MX>(), sp));
    }
  }

  MX MXNode::getRef(const Slice& i, const Slice& j) const {
    return MX::create(new SubRef(shared_from_this<MX>(), i, j));
  }

  MX MXNode::getAssign(const MX& y, const Slice& i, const Slice& j) const {
    return MX::create(new SubAssign(shared_from_this<MX>(), y, i, j));
  }

  MX MXNode::getUnary(int op) const {
    if (operation_checker<F0XChecker>(op) && is_zero()) {
      // If identically zero
      return MX::zeros(sparsity());
    } else {
      // Create a new node
      return MX::create(new UnaryMX(Operation(op), shared_from_this<MX>()));
    }
  }

  MX MXNode::getBinarySwitch(int op, const MX& y) const {
    // Make sure that dimensions match
    casadi_assert_message(sparsity().is_scalar() || y.is_scalar() || sparsity().size()==y.size(),
                          "Dimension mismatch." << "lhs is " << sparsity().dim()
                          << ", while rhs is " << y.dim());

    // Create binary node
    if (sparsity().is_scalar(false)) {
      if (nnz()==0) {
        return toMatrix(MX(0)->getBinary(op, y, true, false), y.sparsity());
      } else {
        return toMatrix(getBinary(op, y, true, false), y.sparsity());
      }
    } else if (y.is_scalar()) {
      if (y.nnz()==0) {
        return toMatrix(getBinary(op, MX(0), false, true), sparsity());
      } else {
        return toMatrix(getBinary(op, y, false, true), sparsity());
      }
    } else {
      casadi_assert_message(sparsity().size() == y.sparsity().size(), "Dimension mismatch.");
      if (sparsity()==y.sparsity()) {
        // Matching sparsities
        return getBinary(op, y, false, false);
      } else {
        // Get the sparsity pattern of the result
        // (ignoring structural zeros giving rise to nonzero result)
        const Sparsity& x_sp = sparsity();
        const Sparsity& y_sp = y.sparsity();
        Sparsity r_sp = x_sp.combine(y_sp, operation_checker<F0XChecker>(op),
                                            operation_checker<FX0Checker>(op));

        // Project the arguments to this sparsity
        MX xx = project(shared_from_this<MX>(), r_sp);
        MX yy = project(y, r_sp);
        return xx->getBinary(op, yy, false, false);
      }
    }
  }

  MX MXNode::getBinary(int op, const MX& y, bool scX, bool scY) const {
    casadi_assert(sparsity()==y.sparsity() || scX || scY);

    if (GlobalOptions::simplification_on_the_fly) {

      // If identically zero due to one argumebt being zero
      if ((operation_checker<F0XChecker>(op) && is_zero()) ||
         (operation_checker<FX0Checker>(op) && y->is_zero())) {
        return MX::zeros(sparsity());
      }

      // Handle special operations (independent of type)
      switch (op) {
      case OP_ADD:
        if (MXNode::is_equal(y.get(), this, maxDepth())) return getUnary(OP_TWICE);
        break;
      case OP_SUB:
      case OP_NE:
      case OP_LT:
        if (MXNode::is_equal(y.get(), this, maxDepth())) return MX::zeros(sparsity());
        break;
      case OP_DIV:
        if (y->is_zero()) return MX::nan(sparsity());
        // fall-through
      case OP_EQ:
      case OP_LE:
        if (MXNode::is_equal(y.get(), this, maxDepth())) return MX::ones(sparsity());
        break;
      case OP_MUL:
        if (MXNode::is_equal(y.get(), this, maxDepth())) return getUnary(OP_SQ);
        break;
      default: break; // no rule
      }

      // Handle special cases for the second argument
      switch (y->op()) {
      case OP_CONST:
        // Make the constant the first argument, if possible
        if (this->op()!=OP_CONST && operation_checker<CommChecker>(op)) {
          return y->getBinary(op, shared_from_this<MX>(), scY, scX);
        } else {
          switch (op) {
          case OP_POW:
            return getBinary(OP_CONSTPOW, y, scX, scY);
          case OP_CONSTPOW:
            if (y->isValue(-1)) return getUnary(OP_INV);
            else if (y->isValue(0)) return MX::ones(sparsity());
            else if (y->isValue(1)) return shared_from_this<MX>();
            else if (y->isValue(2)) return getUnary(OP_SQ);
            break;
          case OP_ADD:
          case OP_SUB:
            if (y->is_zero())
                return scX ? repmat(shared_from_this<MX>(), y.size()) : shared_from_this<MX>();
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
      if (y.is_dense() || operation_checker<FX0Checker>(op)) {
        // Loop over nonzeros
        return MX::create(new BinaryMX<true, false>(Operation(op), shared_from_this<MX>(), y));
      } else {
        // Put a densification node in between
        return getBinary(op, densify(y), true, false);
      }
    } else if (scY) {
      // Check if it is ok to loop over nonzeros only
      if (sparsity().is_dense() || operation_checker<F0XChecker>(op)) {
        // Loop over nonzeros
        return MX::create(new BinaryMX<false, true>(Operation(op), shared_from_this<MX>(), y));
      } else {
        // Put a densification node in between
        return densify(shared_from_this<MX>())->getBinary(op, y, false, true);
      }
    } else {
      // Loop over nonzeros only
      MX rr = MX::create(new BinaryMX<false, false>(Operation(op), shared_from_this<MX>(), y));

      // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
      if (!rr.is_dense() && !operation_checker<F00Checker>(op)) {
        // Get the value for the structural zeros
        double fcn_0(0);
        casadi_math<double>::fun(op, 0, 0, fcn_0);
        rr = densify(rr, fcn_0);
      }
      return rr;
    }
  }

  Matrix<int> MXNode::mapping() const {
    throw CasadiException(string("MXNode::mapping not defined for class ") + typeid(*this).name());
  }

  bool MXNode::sameOpAndDeps(const MXNode* node, int depth) const {
    if (op()!=node->op() || ndep()!=node->ndep())
      return false;
    for (int i=0; i<ndep(); ++i) {
      if (!MX::is_equal(dep(i), node->dep(i), depth-1))
        return false;
    }
    return true;
  }

  MX MXNode::getAssertion(const MX& y, const std::string& fail_message) const {
    return MX::create(new Assertion(shared_from_this<MX>(), y, fail_message));
  }

  MX MXNode::getMonitor(const std::string& comment) const {
    if (nnz()==0) {
      return shared_from_this<MX>();
    } else {
      return MX::create(new Monitor(shared_from_this<MX>(), comment));
    }
  }

  MX MXNode::getFind() const {
    return MX::create(new Find(shared_from_this<MX>()));
  }

  MX MXNode::getDeterminant() const {
    return MX::create(new Determinant(shared_from_this<MX>()));
  }

  MX MXNode::getInverse() const {
    return MX::create(new Inverse(shared_from_this<MX>()));
  }


  MX MXNode::getDot(const MX& y) const {
    casadi_assert_message(
      size2()==y.size2() && size1()==y.size1(),
      "MXNode::dot: Dimension mismatch. dot requires its "
      "two arguments to have equal shapes, but got ("
      << size2() << ", " << size1() << ") and ("
      << y.size2() << ", " << y.size1() << ").");
    if (sparsity()==y.sparsity()) {
      if (sparsity().nnz()==0) {
        return 0;
      } else if (sparsity().is_scalar()) {
        return getBinarySwitch(OP_MUL, y);
      } else {
        return MX::create(new Dot(shared_from_this<MX>(), y));
      }
    } else {
      // Project to pattern intersection
      Sparsity sp = sparsity().intersect(y.sparsity());
      MX xx = project(shared_from_this<MX>(), sp);
      MX yy = project(y, sp);
      return xx->getDot(yy);
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

  MX MXNode::getHorzcat(const vector<MX>& x) const {
    // Check if there is any existing horzcat operation
    for (auto i=x.begin(); i!=x.end(); ++i) {
      if (i->op()==OP_HORZCAT) {
        // Split up
        vector<MX> x_split(x.begin(), i);
        for (; i!=x.end(); ++i) {
          if (i->op()==OP_HORZCAT) {
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

  MX MXNode::get_diagcat(const vector<MX>& x) const {
    // Create a Horzcat node
    return MX::create(new Diagcat(x));
  }

  MX MXNode::getVertcat(const vector<MX>& x) const {
    // Check if there is any existing vertcat operation
    for (auto i=x.begin(); i!=x.end(); ++i) {
      if (i->op()==OP_VERTCAT) {
        // Split up
        vector<MX> x_split(x.begin(), i);
        for (; i!=x.end(); ++i) {
          if (i->op()==OP_VERTCAT) {
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

  vector<MX> MXNode::getHorzsplit(const vector<int>& output_offset) const {
    if (is_zero()) {
      vector<MX> ret =
          MX::createMultipleOutput(new Horzsplit(shared_from_this<MX>(), output_offset));
      for (int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    vector<MX> ret =
        MX::createMultipleOutput(new Horzsplit(shared_from_this<MX>(), output_offset));

    if (GlobalOptions::simplification_on_the_fly) {
      // Simplify horzsplit(horzcat)
      if (op()==OP_HORZCAT) {
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

  MX MXNode::getRepmat(int n, int m) const {
    if (n==1) {
      return MX::create(new HorzRepmat(shared_from_this<MX>(), m));
    } else {
      // Fallback to generic_matrix impl
      return GenericMatrix<MX>::repmat(shared_from_this<MX>(), n, m);
    }
  }

  MX MXNode::getRepsum(int n, int m) const {
    if (n==1) {
      return MX::create(new HorzRepsum(shared_from_this<MX>(), m));
    } else {
      // Fallback to generic_matrix impl
      return GenericMatrix<MX>::repsum(shared_from_this<MX>(), n, m);
    }
  }

  vector<MX> MXNode::get_diagsplit(const vector<int>& offset1,
                                       const vector<int>& offset2) const {
    if (is_zero()) {
      vector<MX> ret =
          MX::createMultipleOutput(new Diagsplit(shared_from_this<MX>(), offset1, offset2));
      for (int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    vector<MX> ret =
        MX::createMultipleOutput(new Diagsplit(shared_from_this<MX>(), offset1, offset2));

    return ret;
  }

  vector<MX> MXNode::getVertsplit(const vector<int>& output_offset) const {
    if (is_zero()) {
      vector<MX> ret =
          MX::createMultipleOutput(new Vertsplit(shared_from_this<MX>(), output_offset));
      for (int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    vector<MX> ret =
        MX::createMultipleOutput(new Vertsplit(shared_from_this<MX>(), output_offset));

    if (GlobalOptions::simplification_on_the_fly) {
      // Simplify vertsplit(vertcat)
      if (op()==OP_VERTCAT) {
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

  void MXNode::copyFwd(const bvec_t* arg, bvec_t* res, int len) {
    if (arg!=res) {
      copy(arg, arg+len, res);
    }
  }

  void MXNode::copyAdj(bvec_t* arg, bvec_t* res, int len) {
    if (arg!=res) {
      for (int k=0; k<len; ++k) {
        *arg++ |= *res;
        *res++ = 0;
      }
    }
  }

  bool MXNode::is_equal(const MXNode* x, const MXNode* y, int depth) {
    if (x==y) {
      return true;
    } else if (depth>0) {
      return x->is_equal(y, depth);
    } else {
      return false;
    }
  }

} // namespace casadi
