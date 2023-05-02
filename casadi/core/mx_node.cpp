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


#include "mx_node.hpp"
#include "casadi_misc.hpp"
#include "transpose.hpp"
#include "reshape.hpp"
#include "sparsity_cast.hpp"
#include "multiplication.hpp"
#include "bilin.hpp"
#include "rank1.hpp"
#include "subref.hpp"
#include "subassign.hpp"
#include "getnonzeros.hpp"
#include "getnonzeros_param.hpp"
#include "setnonzeros.hpp"
#include "setnonzeros_param.hpp"
#include "project.hpp"
#include "solve.hpp"
#include "unary_mx.hpp"
#include "binary_mx.hpp"
#include "determinant.hpp"
#include "inverse.hpp"
#include "dot.hpp"
#include "norm.hpp"
#include "mmin.hpp"
#include "concat.hpp"
#include "split.hpp"
#include "assertion.hpp"
#include "monitor.hpp"
#include "repmat.hpp"
#include "casadi_find.hpp"
#include "casadi_low.hpp"
#include "einstein.hpp"
#include "io_instruction.hpp"
#include "symbolic_mx.hpp"
#include "constant_mx.hpp"
#include "map.hpp"
#include "bspline.hpp"
#include "convexify.hpp"
#include "logsumexp.hpp"

// Template implementations
#include "setnonzeros_impl.hpp"
#include "setnonzeros_param_impl.hpp"
#include "solve_impl.hpp"
#include "binary_mx_impl.hpp"

#include "serializing_stream.hpp"

#include <typeinfo>

namespace casadi {

  MXNode::MXNode() {
    temp = 0;
  }


  MXNode::~MXNode() {

    // Start destruction method if any of the dependencies has dependencies
    for (std::vector<MX>::iterator cc=dep_.begin(); cc!=dep_.end(); ++cc) {
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
          while (!t->dep_.empty()) {
            const MX& ii = t->dep_.back();

            // Skip if constant
            if (ii.is_constant()) {
              t->dep_.pop_back();
              continue;
            }
            // Check if this is the only reference to the element
            if (ii.getCount()==1) {
              // Remove and add to stack
              deletion_stack.push(ii);
              t->dep_.pop_back();
              found_dep = true;
              break;
            } else {
              t->dep_.pop_back();
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

  casadi_int MXNode::n_primitives() const {
    return 1;
  }

  bool MXNode::has_duplicates() const {
    casadi_error("'has_duplicates' not defined for class " + class_name());
  }

  void MXNode::reset_input() const {
    casadi_error("'reset_input' not defined for class " + class_name());
  }

  void MXNode::primitives(std::vector<MX>::iterator& it) const {
    *it++ = shared_from_this<MX>();
  }

  void MXNode::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    *it++ = x;
  }

  MX MXNode::join_primitives(std::vector<MX>::const_iterator& it) const {
    MX ret = *it++;
    if (ret.size()==size()) {
      return ret;
    } else {
      casadi_assert_dev(ret.is_empty(true));
      return MX(size());
    }
  }

  const std::string& MXNode::name() const {
    casadi_error("'name' not defined for class " + class_name());
  }

  std::string MXNode::class_name() const {
    // Lazy solution
    return typeid(*this).name();
  }

  bool MXNode::__nonzero__() const {
    casadi_error("Can only determine truth value of a numeric MX.");
  }

  casadi_int MXNode::n_dep() const {
    return dep_.size();
  }

  casadi_int MXNode::ind() const {
    casadi_error("'ind' not defined for class " + class_name());
  }

  casadi_int MXNode::segment() const {
    casadi_error("'segment' not defined for class " + class_name());
  }

  casadi_int MXNode::offset() const {
    casadi_error("'offset' not defined for class " + class_name());
  }

  void MXNode::set_sparsity(const Sparsity& sparsity) {
    sparsity_ = sparsity;
  }

  void MXNode::set_dep(const MX& dep) {
    dep_.resize(1);
    dep_[0] = dep;
  }

  void MXNode::set_dep(const MX& dep1, const MX& dep2) {
    dep_.resize(2);
    dep_[0] = dep1;
    dep_[1] = dep2;
  }

  void MXNode::set_dep(const MX& dep1, const MX& dep2, const MX& dep3) {
    dep_.resize(3);
    dep_[0] = dep1;
    dep_[1] = dep2;
    dep_[2] = dep3;
  }

  void MXNode::set_dep(const std::vector<MX>& dep) {
    dep_ = dep;
  }

  const Sparsity& MXNode::sparsity(casadi_int oind) const {
    casadi_assert(oind==0, "Index out of bounds");
    return sparsity_;
  }

  void MXNode::disp(std::ostream& stream, bool more) const {
    // Find out which noded can be inlined
    std::map<const MXNode*, casadi_int> nodeind;
    can_inline(nodeind);

    // Print expression
    std::vector<std::string> intermed;
    std::string s = print_compact(nodeind, intermed);

    // Print intermediate expressions
    for (casadi_int i=0; i<intermed.size(); ++i)
      stream << "@" << (i+1) << "=" << intermed[i] << ", ";

    // Print this
    stream << s;
  }

  void MXNode::can_inline(std::map<const MXNode*, casadi_int>& nodeind) const {
    // Add or mark node in map
    std::map<const MXNode*, casadi_int>::iterator it=nodeind.find(this);
    if (it==nodeind.end()) {
      // First time encountered, mark inlined
      nodeind.insert(it, std::make_pair(this, 0));

      // Handle dependencies with recursion
      for (casadi_int i=0; i<n_dep(); ++i) {
        dep(i)->can_inline(nodeind);
      }
    } else if (it->second==0 && op()!=OP_PARAMETER) {
      // Node encountered before, do not inline (except if symbolic primitive)
      it->second = -1;
    }
  }

  std::string MXNode::print_compact(std::map<const MXNode*, casadi_int>& nodeind,
                                   std::vector<std::string>& intermed) const {
    // Get reference to node index
    casadi_int& ind = nodeind[this];

    // If positive, already in intermediate expressions
    if (ind>0) return "@" + str(ind);

    // Get expressions for dependencies
    std::vector<std::string> arg(n_dep());
    for (casadi_int i=0; i<arg.size(); ++i) {
      arg[i] = dep(i)->print_compact(nodeind, intermed);
    }

    // Get expression for this
    std::string s = disp(arg);

    // Decide what to do with the expression
    if (ind==0) {
      // Inline expression
      return s;
    } else {
      // Add to list of intermediate expressions and return reference
      intermed.push_back(s);
      ind = intermed.size(); // For subsequent references
      return "@" + str(ind);
    }
  }

  const Function& MXNode::which_function() const {
    casadi_error("'which_function' not defined for class " + class_name());
  }

  casadi_int MXNode::which_output() const {
    casadi_error("'which_output' not defined for class " + class_name());
  }

  int MXNode::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    casadi_error("'eval' not defined for class " + class_name());
    return 1;
  }

  int MXNode::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    casadi_error("'eval_sx' not defined for class " + class_name());
    return 1;
  }

  void MXNode::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    casadi_error("'eval_mx' not defined for class " + class_name());
  }

  void MXNode::ad_forward(const std::vector<std::vector<MX> >& fseed,
                       std::vector<std::vector<MX> >& fsens) const {
    casadi_error("'ad_forward' not defined for class " + class_name());
  }

  void MXNode::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                       std::vector<std::vector<MX> >& asens) const {
    casadi_error("'ad_reverse' not defined for class " + class_name());
  }

  int MXNode::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // By default, everything depends on everything
    bvec_t all_depend(0);

    // Get dependencies of all inputs
    for (casadi_int k=0; k<n_dep(); ++k) {
      const bvec_t* v = arg[k];
      for (casadi_int i=0; i<dep(k).nnz(); ++i) {
        all_depend |= v[i];
      }
    }

    // Propagate to all outputs
    for (casadi_int k=0; k<nout(); ++k) {
      bvec_t* v = res[k];
      for (casadi_int i=0; i<sparsity(k).nnz(); ++i) {
        v[i] = all_depend;
      }
    }
    return 0;
  }

  int MXNode::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // By default, everything depends on everything
    bvec_t all_depend(0);

    // Get dependencies of all outputs
    for (casadi_int k=0; k<nout(); ++k) {
      bvec_t* v = res[k];
      for (casadi_int i=0; i<sparsity(k).nnz(); ++i) {
        all_depend |= v[i];
        v[i] = 0;
      }
    }

    // Propagate to all inputs
    for (casadi_int k=0; k<n_dep(); ++k) {
      bvec_t* v = arg[k];
      for (casadi_int i=0; i<dep(k).nnz(); ++i) {
        v[i] |= all_depend;
      }
    }
    return 0;
  }

  MX MXNode::get_output(casadi_int oind) const {
    casadi_assert(oind==0, "Output index out of bounds");
    return shared_from_this<MX>();
  }

  void MXNode::generate(CodeGenerator& g,
      const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    casadi_warning("Cannot code generate MX nodes of type " + class_name() +
                   "The generation will proceed, but compilation of the code will "
                   "not be possible.");
    g << "#error " <<  class_name() << ": " << arg << " => " << res << '\n';
  }

  double MXNode::to_double() const {
    casadi_error("'to_double' not defined for class " + class_name());
  }

  DM MXNode::get_DM() const {
    casadi_error("'get_DM' not defined for class " + class_name());
  }

  MX MXNode::get_transpose() const {
    if (sparsity().is_scalar()) {
      return shared_from_this<MX>();
    } else if (sparsity().is_vector()) {
      return get_reshape(sparsity().T());
    } else if (sparsity().is_dense()) {
      return MX::create(new DenseTranspose(shared_from_this<MX>()));
    } else {
      return MX::create(new Transpose(shared_from_this<MX>()));
    }
  }

  MX MXNode::get_reshape(const Sparsity& sp) const {
    casadi_assert_dev(sp.is_reshape(sparsity()));
    if (sp==sparsity()) {
      return shared_from_this<MX>();
    } else {
      return MX::create(new Reshape(shared_from_this<MX>(), sp));
    }
  }

  MX MXNode::get_sparsity_cast(const Sparsity& sp) const {
    casadi_assert_dev(sp.nnz()==nnz());
    if (sp==sparsity()) {
      return shared_from_this<MX>();
    } else {
      return MX::create(new SparsityCast(shared_from_this<MX>(), sp));
    }
  }

  Dict MXNode::info() const {
    return Dict();
  }

  void MXNode::serialize(SerializingStream& s) const {
    serialize_type(s);
    serialize_body(s);
  }

  void MXNode::serialize_body(SerializingStream& s) const {
    s.pack("MXNode::deps", dep_);
    s.pack("MXNode::sp", sparsity_);
  }

  void MXNode::serialize_type(SerializingStream& s) const {
    s.pack("MXNode::op", static_cast<int>(op()));
  }

  MXNode::MXNode(DeserializingStream& s) {
    temp = 0;

    s.unpack("MXNode::deps", dep_);
    s.unpack("MXNode::sp", sparsity_);
  }


  MXNode* MXNode::deserialize(DeserializingStream& s) {
    int op;
    s.unpack("MXNode::op", op);

    if (casadi_math<MX>::is_binary(op)) {
      return BinaryMX<false, false>::deserialize(s);
    } else if (casadi_math<MX>::is_unary(op)) {
      return UnaryMX::deserialize(s);
    }

    auto it = MXNode::deserialize_map.find(op);
    if (it==MXNode::deserialize_map.end()) {
      casadi_error("Not implemented op " + str(casadi_int(op)) + ":" + str(OP_GETNONZEROS));
    } else {
      return it->second(s);
    }
  }


  MX MXNode::get_mac(const MX& y, const MX& z) const {
    if (sparsity().is_orthonormal() && y.is_column() && y.is_dense()
        && y.sparsity()==z.sparsity() && z.is_zero()) {
      std::vector<casadi_int> perm = sparsity().permutation_vector();
      MX nz = sparsity_cast(shared_from_this<MX>(), Sparsity::dense(nnz()));
      return (nz*y)(perm);
    }
    // Get reference to transposed first argument
    MX x = shared_from_this<MX>();

    casadi_assert(y.size2()==z.size2(),
      "Dimension error x.mac(z). Got y=" + str(y.size2()) + " and z=" + z.dim() + ".");
    casadi_assert(x.size1()==z.size1(),
      "Dimension error x.mac(z). Got x=" + x.dim() + " and z=" + z.dim() + ".");
    casadi_assert(y.size1()==x.size2(),
      "Dimension error x.mac(z). Got y=" + str(y.size1()) + " and x" + x.dim() + ".");
    if (x.is_dense() && y.is_dense() && z.is_dense()) {
      return MX::create(new DenseMultiplication(z, x, y));
    } else {
      return MX::create(new Multiplication(z, x, y));
    }
  }

  MX MXNode::get_einstein(const MX& A, const MX& B,
      const std::vector<casadi_int>& dim_c, const std::vector<casadi_int>& dim_a,
      const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& c, const std::vector<casadi_int>& a,
      const std::vector<casadi_int>& b) const {

    if (A.is_zero() || B.is_zero())
      return shared_from_this<MX>();

    MX C = densify(shared_from_this<MX>());

    if (A.is_constant() && B.is_constant() && C.is_constant()) {
      // Constant folding
      DM Ac = A->get_DM();
      DM Bc = B->get_DM();
      DM Cc = C->get_DM();
      return einstein(vec(densify(Ac)), vec(densify(Bc)), vec(densify(Cc)),
        dim_a, dim_b, dim_c, a, b, c);
    }

    return MX::create(new Einstein(C, densify(A), densify(B), dim_c, dim_a, dim_b, c, a, b));
  }

  MX MXNode::get_bilin(const MX& x, const MX& y) const {
    return MX::create(new Bilin(shared_from_this<MX>(), x, y));
  }

  MX MXNode::get_rank1(const MX& alpha, const MX& x, const MX& y) const {
    return MX::create(new Rank1(shared_from_this<MX>(), alpha, x, y));
  }

  MX MXNode::get_logsumexp() const {
    return MX::create(new LogSumExp(shared_from_this<MX>()));
  }

  MX MXNode::get_solve_triu(const MX& r, bool tr) const {
    if (tr) {
      return MX::create(new TriuSolve<true>(densify(r), shared_from_this<MX>()));
    } else {
      return MX::create(new TriuSolve<false>(densify(r), shared_from_this<MX>()));
    }
  }

  MX MXNode::get_solve_tril(const MX& r, bool tr) const {
    if (tr) {
      return MX::create(new TrilSolve<true>(densify(r), shared_from_this<MX>()));
    } else {
      return MX::create(new TrilSolve<false>(densify(r), shared_from_this<MX>()));
    }
  }

  MX MXNode::get_solve_triu_unity(const MX& r, bool tr) const {
    if (tr) {
      return MX::create(new TriuSolveUnity<true>(densify(r), shared_from_this<MX>()));
    } else {
      return MX::create(new TriuSolveUnity<false>(densify(r), shared_from_this<MX>()));
    }
  }

  MX MXNode::get_solve_tril_unity(const MX& r, bool tr) const {
    if (tr) {
      return MX::create(new TrilSolveUnity<true>(densify(r), shared_from_this<MX>()));
    } else {
      return MX::create(new TrilSolveUnity<false>(densify(r), shared_from_this<MX>()));
    }
  }

  MX MXNode::get_solve(const MX& r, bool tr, const Linsol& linear_solver) const {
    if (tr) {
      return MX::create(new LinsolCall<true>(densify(r), shared_from_this<MX>(), linear_solver));
    } else {
      return MX::create(new LinsolCall<false>(densify(r), shared_from_this<MX>(), linear_solver));
    }
  }

  MX MXNode::get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const {
    if (sparsity().is_dense() && is_range(nz, 0, nnz())) {
      return sparsity_cast(shared_from_this<MX>(), sp);
    }
    return GetNonzeros::create(sp, shared_from_this<MX>(), nz);
  }

  MX MXNode::get_nz_ref(const MX& nz) const {
    return GetNonzerosParam::create(shared_from_this<MX>(), nz);
  }

  MX MXNode::get_nz_ref(const MX& inner, const Slice& outer) const {
    if (outer.all()==std::vector<casadi_int>{0}) {
      return get_nz_ref(inner);
    } else {
      return GetNonzerosParam::create(shared_from_this<MX>(), inner, outer);
    }
  }

  MX MXNode::get_nz_ref(const Slice& inner, const MX& outer) const {
    if (inner.all()==std::vector<casadi_int>{0}) {
      return get_nz_ref(outer);
    } else {
      return GetNonzerosParam::create(shared_from_this<MX>(), inner, outer);
    }
  }

  MX MXNode::get_nz_ref(const MX& inner, const MX& outer) const {
    return GetNonzerosParam::create(shared_from_this<MX>(), inner, outer);
  }

  MX MXNode::get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const {
    // Check if any element needs to be set at all
    bool set_any = false;
    for (auto i=nz.begin(); i!=nz.end() && !set_any; ++i) {
      set_any = *i >= 0;
    }
    if (!set_any) return y;

    return SetNonzeros<false>::create(y, shared_from_this<MX>(), nz);
  }


  MX MXNode::get_nzadd(const MX& y, const std::vector<casadi_int>& nz) const {
    if (nz.empty() || is_zero()) {
      return y;
    } else {
      return SetNonzeros<true>::create(y, shared_from_this<MX>(), nz);
    }
  }

  MX MXNode::get_nzassign(const MX& y, const MX& nz) const {
    return SetNonzerosParam<false>::create(y, shared_from_this<MX>(), nz);
  }

  MX MXNode::get_nzassign(const MX& y, const MX& inner, const Slice& outer) const {
    return SetNonzerosParam<false>::create(y, shared_from_this<MX>(), inner, outer);
  }

  MX MXNode::get_nzassign(const MX& y, const Slice& inner, const MX& outer) const {
    return SetNonzerosParam<false>::create(y, shared_from_this<MX>(), inner, outer);
  }

  MX MXNode::get_nzassign(const MX& y, const MX& inner, const MX& outer) const {
    return SetNonzerosParam<false>::create(y, shared_from_this<MX>(), inner, outer);
  }

  MX MXNode::get_nzadd(const MX& y, const MX& nz) const {
    if (nz.is_empty() || is_zero()) {
      return y;
    } else {
      return SetNonzerosParam<true>::create(y, shared_from_this<MX>(), nz);
    }
  }

  MX MXNode::get_nzadd(const MX& y, const MX& inner, const Slice& outer) const {
    if (inner.is_empty() || outer.is_empty() || is_zero()) {
      return y;
    } else {
      return SetNonzerosParam<true>::create(y, shared_from_this<MX>(), inner, outer);
    }
  }

  MX MXNode::get_nzadd(const MX& y, const Slice& inner, const MX& outer) const {
    if (outer.is_empty() || outer.is_empty() || is_zero()) {
      return y;
    } else {
      return SetNonzerosParam<true>::create(y, shared_from_this<MX>(), inner, outer);
    }
  }

  MX MXNode::get_nzadd(const MX& y, const MX& inner, const MX& outer) const {
    if (inner.is_empty() || outer.is_empty() || is_zero()) {
      return y;
    } else {
      return SetNonzerosParam<true>::create(y, shared_from_this<MX>(), inner, outer);
    }
  }

  MX MXNode::get_project(const Sparsity& sp) const {
    if (sp==sparsity()) {
      return shared_from_this<MX>();
    } else if (sp.nnz()==0) {
      return MX::zeros(sp);
    } else if (sp.is_dense()) {
      return MX::create(new Densify(shared_from_this<MX>(), sp));
    } else if (sparsity().is_dense()) {
      return MX::create(new Sparsify(shared_from_this<MX>(), sp));
    } else {
      return MX::create(new Project(shared_from_this<MX>(), sp));
    }
  }

  MX MXNode::get_subref(const Slice& i, const Slice& j) const {
    return MX::create(new SubRef(shared_from_this<MX>(), i, j));
  }

  MX MXNode::get_subassign(const MX& y, const Slice& i, const Slice& j) const {
    return MX::create(new SubAssign(shared_from_this<MX>(), y, i, j));
  }

  MX MXNode::get_unary(casadi_int op) const {
    if (operation_checker<F0XChecker>(op) && is_zero()) {
      // If identically zero
      return MX::zeros(sparsity());
    } else {
      // Create a new node
      return MX::create(new UnaryMX(Operation(op), shared_from_this<MX>()));
    }
  }

  MX MXNode::get_binary(casadi_int op, const MX& y) const {
    // If-else-zero nodes are always simplified at top level to avoid NaN propagation
    if (y.op() == OP_IF_ELSE_ZERO) {
      if (op == OP_MUL) {
        // (Rule 1.) x * if_else_zero(c, y), simplified to if_else_zero(c, x * y)
        // Background: x is often a partial derivative and may evaluate to INF or NAN.
        // The simplification ensures that the zero seed corresponding to an inactive branch does
        // not give rise to any NaN contribution to the derivative due to NaN * 0 == NaN.
        return if_else_zero(y.dep(0), shared_from_this<MX>() * y.dep(1));
      } else if (op == OP_ADD && this->op() == OP_IF_ELSE_ZERO && MX::is_equal(dep(0), y.dep(0))) {
        // (Rule 2.) if_else_zero(c, x) + if_else_zero(c, y) is simplified to if_else_zero(c, x + y)
        // Background: During the backward propagation, seeds are added together. Without this rule,
        // the addition node can prevent rule (1.) from working in subsequent steps.
        return if_else_zero(y.dep(0), dep(1) + y.dep(1));
      }
    } else if (this->op() == OP_IF_ELSE_ZERO && op == OP_MUL) {
      // Same as Rule 1. above, but with factors swapped. For symmetry.
      return if_else_zero(dep(0), dep(1) * y);
    }
    // Create binary node
    if (sparsity().is_scalar(false)) {
      if (nnz()==0) {
        if (operation_checker<F0XChecker>(op)) return MX::zeros(Sparsity(y.size()));
        return to_matrix(MX(0)->_get_binary(op, y, true, false), y.sparsity());
      } else {
        return to_matrix(_get_binary(op, y, true, false), y.sparsity());
      }
    } else if (y.is_scalar()) {
      if (y.nnz()==0) {
        if (operation_checker<FX0Checker>(op)) return MX::zeros(Sparsity(size()));
        return to_matrix(_get_binary(op, MX(0), false, true), sparsity());
      } else {
        return to_matrix(_get_binary(op, y, false, true), sparsity());
      }
    } else {
      casadi_assert(sparsity().size() == y.sparsity().size(), "Dimension mismatch.");
      if (sparsity()==y.sparsity()) {
        // Matching sparsities
        return _get_binary(op, y, false, false);
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
        return xx->_get_binary(op, yy, false, false);
      }
    }
  }

  MX MXNode::_get_binary(casadi_int op, const MX& y, bool scX, bool scY) const {
    casadi_assert_dev(sparsity()==y.sparsity() || scX || scY);

    if (GlobalOptions::simplification_on_the_fly) {

      // If identically zero due to one argument being zero
      if ((operation_checker<F0XChecker>(op) && nnz()==0) ||
         (operation_checker<FX0Checker>(op) && y->nnz()==0)) {
        return MX::zeros(Sparsity(size()));
      }

      if ((operation_checker<F0XChecker>(op) && is_zero()) ||
         (operation_checker<FX0Checker>(op) && y->is_zero())) {
        return MX::zeros(sparsity());
      }

      // Handle special operations (independent of type)
      switch (op) {
      case OP_ADD:
        if (MXNode::is_equal(y.get(), this, maxDepth())) return get_unary(OP_TWICE);
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
        if (MXNode::is_equal(y.get(), this, maxDepth())) return get_unary(OP_SQ);
        break;
      default: break; // no rule
      }

      // Handle special cases for the second argument
      switch (y->op()) {
      case OP_CONST:
        // Make the constant the first argument, if possible
        if (this->op()!=OP_CONST && operation_checker<CommChecker>(op)) {
          return y->_get_binary(op, shared_from_this<MX>(), scY, scX);
        } else {
          switch (op) {
          case OP_POW:
            return _get_binary(OP_CONSTPOW, y, scX, scY);
          case OP_CONSTPOW:
            if (y->is_value(-1)) return get_unary(OP_INV);
            else if (y->is_value(0)) return MX::ones(sparsity());
            else if (y->is_value(1)) return shared_from_this<MX>();
            else if (y->is_value(2)) return get_unary(OP_SQ);
            break;
          case OP_ADD:
          case OP_SUB:
            if (y->is_zero())
                return scX ? repmat(shared_from_this<MX>(), y.size()) : shared_from_this<MX>();
            break;
          case OP_MUL:
            if (y->is_value(1)) return shared_from_this<MX>();
            break;
          case OP_DIV:
            if (y->is_value(1)) return shared_from_this<MX>();
            else if (y->is_value(0.5)) return get_unary(OP_TWICE);
            break;
          default: break; // no rule
          }
        }
        break;
      case OP_NEG:
        if (op==OP_ADD) {
          return _get_binary(OP_SUB, y->dep(), scX, scY);
        } else if (op==OP_SUB) {
          return _get_binary(OP_ADD, y->dep(), scX, scY);
        } else if (op==OP_MUL) {
          return -_get_binary(OP_MUL, y->dep(), scX, scY);
        } else if (op==OP_DIV) {
          return -_get_binary(OP_DIV, y->dep(), scX, scY);
        }
        break;
      case OP_INV:
        if (op==OP_MUL) {
          return _get_binary(OP_DIV, y->dep(), scX, scY);
        } else if (op==OP_DIV) {
          return _get_binary(OP_MUL, y->dep(), scX, scY);
        }
        break;
      default: break; // no rule
      }

    }

    if (scX) {
      // Check if it is ok to loop over nonzeros only
      if (y.is_dense() || operation_checker<FX0Checker>(op) ||
          (is_zero() && operation_checker<F00Checker>(op))) {
        // Loop over nonzeros
        return MX::create(new BinaryMX<true, false>(Operation(op), shared_from_this<MX>(), y));
      } else {
        // Put a densification node in between
        return _get_binary(op, densify(y), true, false);
      }
    } else if (scY) {
      // Check if it is ok to loop over nonzeros only
      if (sparsity().is_dense() || operation_checker<F0XChecker>(op) ||
          (y.is_zero() && operation_checker<F00Checker>(op))) {
        // Loop over nonzeros
        return MX::create(new BinaryMX<false, true>(Operation(op), shared_from_this<MX>(), y));
      } else {
        // Put a densification node in between
        return densify(shared_from_this<MX>())->_get_binary(op, y, false, true);
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

  Matrix<casadi_int> MXNode::mapping() const {
    casadi_error("'mapping' not defined for class " + class_name());
  }

  bool MXNode::sameOpAndDeps(const MXNode* node, casadi_int depth) const {
    if (op()!=node->op() || n_dep()!=node->n_dep())
      return false;
    for (casadi_int i=0; i<n_dep(); ++i) {
      if (!MX::is_equal(dep(i), node->dep(i), depth-1))
        return false;
    }
    return true;
  }

  MX MXNode::get_assert(const MX& y, const std::string& fail_message) const {
    return MX::create(new Assertion(shared_from_this<MX>(), y, fail_message));
  }

  MX MXNode::get_monitor(const std::string& comment) const {
    if (nnz()==0) {
      return shared_from_this<MX>();
    } else {
      return MX::create(new Monitor(shared_from_this<MX>(), comment));
    }
  }

  MX MXNode::get_find() const {
    MX x = shared_from_this<MX>();
    casadi_assert(x.is_vector(), "Argument must be vector, got " + x.dim() + ".");
    if (x.is_column()) {
      return MX::create(new Find(shared_from_this<MX>()));
    } else {
      return find(x.T());
    }
  }

  MX MXNode::get_low(const MX& v, const Dict& options) const {
    return MX::create(new Low(v, shared_from_this<MX>(), options));
  }

  MX MXNode::get_bspline(const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<double>& coeffs,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode) const {
    MX x = shared_from_this<MX>();
    return MX::create(new BSpline(x, knots, offset, coeffs, degree, m, lookup_mode));
  }

  MX MXNode::get_bspline(const MX& coeffs,
            const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode) const {
    MX x = shared_from_this<MX>();
    return MX::create(new BSplineParametric(x, coeffs, knots, offset, degree, m, lookup_mode));
  }

  MX MXNode::get_convexify(const Dict& opts) const {
    return MX::create(new Convexify(shared_from_this<MX>(), opts));
  }

  MX MXNode::get_det() const {
    return MX::create(new Determinant(shared_from_this<MX>()));
  }

  MX MXNode::get_inv() const {
    return MX::create(new Inverse(shared_from_this<MX>()));
  }


  MX MXNode::get_dot(const MX& y) const {
    casadi_assert(
      size2()==y.size2() && size1()==y.size1(),
      "MXNode::dot: Dimension mismatch. dot requires its "
      "two arguments to have equal shapes, but got ("
      + str(size2()) + ", " + str(size1()) + ") and ("
      + str(y.size2()) + ", " + str(y.size1()) + ").");
    if (sparsity()==y.sparsity()) {
      if (sparsity().nnz()==0) {
        return 0;
      } else if (sparsity().is_scalar()) {
        return get_binary(OP_MUL, y);
      } else {
        return MX::create(new Dot(shared_from_this<MX>(), y));
      }
    } else {
      // Project to pattern intersection
      Sparsity sp = sparsity().intersect(y.sparsity());
      MX xx = project(shared_from_this<MX>(), sp);
      MX yy = project(y, sp);
      return xx->get_dot(yy);
    }
  }

  MX MXNode::get_norm_fro() const {
    return MX::create(new NormF(shared_from_this<MX>()));
  }

  MX MXNode::get_norm_2() const {
    return MX::create(new Norm2(shared_from_this<MX>()));
  }

  MX MXNode::get_norm_inf() const {
    return MX::create(new NormInf(shared_from_this<MX>()));
  }

  MX MXNode::get_norm_1() const {
    return MX::create(new Norm1(shared_from_this<MX>()));
  }

  MX MXNode::get_mmin() const {
    if (sparsity_.is_empty()) return MX();
    return MX::create(new MMin(shared_from_this<MX>()));
  }

  MX MXNode::get_mmax() const {
    if (sparsity_.is_empty()) return MX();
    return MX::create(new MMax(shared_from_this<MX>()));
  }

  MX MXNode::get_horzcat(const std::vector<MX>& x) const {
    // Check if there is any existing horzcat operation
    for (auto i=x.begin(); i!=x.end(); ++i) {
      if (i->op()==OP_HORZCAT) {
        // Split up
        std::vector<MX> x_split(x.begin(), i);
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

  MX MXNode::get_diagcat(const std::vector<MX>& x) const {
    // Create a Horzcat node
    return MX::create(new Diagcat(x));
  }

  MX MXNode::get_vertcat(const std::vector<MX>& x) const {
    // Check if there is any existing vertcat operation
    for (auto i=x.begin(); i!=x.end(); ++i) {
      if (i->op()==OP_VERTCAT) {
        // Split up
        std::vector<MX> x_split(x.begin(), i);
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

  std::vector<MX> MXNode::get_horzsplit(const std::vector<casadi_int>& output_offset) const {
    if (is_zero()) {
      std::vector<MX> ret =
          MX::createMultipleOutput(new Horzsplit(shared_from_this<MX>(), output_offset));
      for (casadi_int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    std::vector<MX> ret =
        MX::createMultipleOutput(new Horzsplit(shared_from_this<MX>(), output_offset));

    if (GlobalOptions::simplification_on_the_fly) {
      // Simplify horzsplit(horzcat)
      if (op()==OP_HORZCAT) {
        casadi_int offset_deps = 0;
        casadi_int j = 0;
        for (casadi_int i=0;i<output_offset.size();++i) {
          while (offset_deps<output_offset[i]) { offset_deps+=dep(j).size2();++j; }
          if (j>=n_dep()) j = n_dep()-1;
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

  MX MXNode::get_repmat(casadi_int n, casadi_int m) const {
    if (n==1) {
      return MX::create(new HorzRepmat(shared_from_this<MX>(), m));
    } else {
      // Fallback to generic_matrix impl
      return GenericMatrix<MX>::repmat(shared_from_this<MX>(), n, m);
    }
  }

  MX MXNode::get_repsum(casadi_int n, casadi_int m) const {
    if (n==1) {
      return MX::create(new HorzRepsum(shared_from_this<MX>(), m));
    } else {
      // Fallback to generic_matrix impl
      return GenericMatrix<MX>::repsum(shared_from_this<MX>(), n, m);
    }
  }

  std::vector<MX> MXNode::get_diagsplit(const std::vector<casadi_int>& offset1,
                                       const std::vector<casadi_int>& offset2) const {
    if (is_zero()) {
      std::vector<MX> ret =
          MX::createMultipleOutput(new Diagsplit(shared_from_this<MX>(), offset1, offset2));
      for (casadi_int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    std::vector<MX> ret =
        MX::createMultipleOutput(new Diagsplit(shared_from_this<MX>(), offset1, offset2));

    return ret;
  }

  std::vector<MX> MXNode::get_vertsplit(const std::vector<casadi_int>& output_offset) const {
    if (is_zero()) {
      std::vector<MX> ret =
          MX::createMultipleOutput(new Vertsplit(shared_from_this<MX>(), output_offset));
      for (casadi_int i=0;i<ret.size();++i) {
        ret[i]=MX::zeros(ret[i].sparsity());
      }
      return ret;
    }
    std::vector<MX> ret =
        MX::createMultipleOutput(new Vertsplit(shared_from_this<MX>(), output_offset));

    if (GlobalOptions::simplification_on_the_fly) {
      // Simplify vertsplit(vertcat)
      if (op()==OP_VERTCAT) {
        casadi_int offset_deps = 0;
        casadi_int j = 0;
        for (casadi_int i=0;i<output_offset.size();++i) {
          while (offset_deps<output_offset[i]) { offset_deps+=dep(j).size1();++j; }
          if (j>=n_dep()) j = n_dep()-1;
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

  void MXNode::copy_fwd(const bvec_t* arg, bvec_t* res, casadi_int len) {
    if (arg!=res) {
      std::copy(arg, arg+len, res);
    }
  }

  void MXNode::copy_rev(bvec_t* arg, bvec_t* res, casadi_int len) {
    if (arg!=res) {
      for (casadi_int k=0; k<len; ++k) {
        *arg++ |= *res;
        *res++ = 0;
      }
    }
  }

  bool MXNode::is_equal(const MXNode* x, const MXNode* y, casadi_int depth) {
    if (x==y) {
      return true;
    } else if (depth>0) {
      return x->is_equal(y, depth);
    } else {
      return false;
    }
  }


  // Note: binary/unary operations are omitted here
  std::map<casadi_int, MXNode* (*)(DeserializingStream&)> MXNode::deserialize_map = {
    {OP_INPUT, Input::deserialize},
    {OP_OUTPUT, Output::deserialize},
    {OP_PARAMETER, SymbolicMX::deserialize},
    {OP_CONST, ConstantMX::deserialize},
    {OP_CALL, Call::deserialize},
    {OP_FIND, Find::deserialize},
    {OP_LOW, Low::deserialize},
    //{OP_MAP, Map::deserialize}, Map is a function
    {OP_MTIMES, Multiplication::deserialize},
    {OP_SOLVE, LinsolCall<false>::deserialize},
    {OP_TRANSPOSE, Transpose::deserialize},
    {OP_DETERMINANT, Determinant::deserialize},
    {OP_INVERSE, Inverse::deserialize},
    {OP_DOT, Dot::deserialize},
    {OP_BILIN, Bilin::deserialize},
    {OP_RANK1, Rank1::deserialize},
    {OP_HORZCAT, Horzcat::deserialize},
    {OP_VERTCAT, Vertcat::deserialize},
    {OP_DIAGCAT, Diagcat::deserialize},
    {OP_HORZSPLIT, Horzsplit::deserialize},
    {OP_VERTSPLIT, Vertsplit::deserialize},
    {OP_DIAGSPLIT, Diagsplit::deserialize},
    {OP_RESHAPE, Reshape::deserialize},
    {OP_SPARSITY_CAST, SparsityCast::deserialize},
    // OP_SUBREF
    // OP_SUBASSIGN,
    {OP_GETNONZEROS, GetNonzeros::deserialize},
    {OP_GETNONZEROS_PARAM, GetNonzerosParam::deserialize},
    {OP_ADDNONZEROS, SetNonzeros<true>::deserialize},
    {OP_ADDNONZEROS_PARAM, SetNonzerosParam<true>::deserialize},
    {OP_SETNONZEROS, SetNonzeros<false>::deserialize},
    {OP_SETNONZEROS_PARAM, SetNonzerosParam<false>::deserialize},
    {OP_PROJECT, Project::deserialize},
    {OP_ASSERTION, Assertion::deserialize},
    {OP_MONITOR, Monitor::deserialize},
    {OP_NORM1, Norm1::deserialize},
    {OP_NORM2, Norm2::deserialize},
    {OP_NORMINF, NormInf::deserialize},
    {OP_NORMF, NormF::deserialize},
    {OP_MMIN, MMin::deserialize},
    {OP_MMAX, MMax::deserialize},
    {OP_HORZREPMAT, HorzRepmat::deserialize},
    {OP_HORZREPSUM, HorzRepsum::deserialize},
    //OP_ERFINV,
    //OP_PRINTME,
    //OP_LIFT,
    //OP_EINSTEIN
    {OP_BSPLINE, BSplineCommon::deserialize},
    {OP_CONVEXIFY, Convexify::deserialize},
    {OP_LOGSUMEXP, LogSumExp::deserialize},
    {-1, OutputNode::deserialize}
  };


} // namespace casadi
