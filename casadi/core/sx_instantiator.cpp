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

#define CASADI_SX_INSTANTIATOR_CPP
#include "matrix_impl.hpp"

#include "sx_function.hpp"

using namespace std;

namespace casadi {

  template<>
  bool CASADI_EXPORT SX::__nonzero__() const {
    casadi_assert(numel()==1,
      "Only scalar Matrix could have a truth value, but you "
      "provided a shape" + dim());
    return nonzeros().at(0).__nonzero__();
  }

  template<>
  void CASADI_EXPORT SX::set_max_depth(casadi_int eq_depth) {
    SXNode::eq_depth_ = eq_depth;
  }

  template<>
  casadi_int CASADI_EXPORT SX::get_max_depth() {
    return SXNode::eq_depth_;
  }

  template<>
  SX CASADI_EXPORT SX::_sym(const string& name, const Sparsity& sp) {
    // Create a dense n-by-m matrix
    vector<SXElem> retv;

    // Check if individial names have been provided
    if (name[0]=='[') {

      // Make a copy of the string and modify it as to remove the special characters
      string modname = name;
      for (string::iterator it=modname.begin(); it!=modname.end(); ++it) {
        switch (*it) {
        case '(': case ')': case '[': case ']': case '{': case '}': case ',': case ';': *it = ' ';
        }
      }

      istringstream iss(modname);
      string varname;

      // Loop over elements
      while (!iss.fail()) {
        // Read the name
        iss >> varname;

        // Append to the return vector
        if (!iss.fail())
          retv.push_back(SXElem::sym(varname));
      }
    } else if (sp.is_scalar(true)) {
      retv.push_back(SXElem::sym(name));
    } else {
      // Scalar
      stringstream ss;
      for (casadi_int k=0; k<sp.nnz(); ++k) {
        ss.str("");
        ss << name << "_" << k;
        retv.push_back(SXElem::sym(ss.str()));
      }
    }

    // Determine dimensions automatically if empty
    if (sp.is_scalar(true)) {
      return SX(retv);
    } else {
      return SX(sp, retv, false);
    }
  }

  template<>
  bool CASADI_EXPORT SX::is_regular() const {
    // First pass: ignore symbolics
    for (casadi_int i=0; i<nnz(); ++i) {
      const SXElem& x = nonzeros().at(i);
      if (x.is_constant()) {
        if (x.is_nan() || x.is_inf() || x.is_minus_inf()) return false;
      }
    }
    // Second pass: don't ignore symbolics
    for (casadi_int i=0; i<nnz(); ++i) {
      if (!nonzeros().at(i).is_regular()) return false;
    }
    return true;
  }

  template<>
  bool CASADI_EXPORT SX::is_smooth() const {
    // Make a function
    Function temp("temp", {SX()}, {*this});

    // Run the function on the temporary variable
    SXFunction* t = temp.get<SXFunction>();
    return t->is_smooth();
  }

  template<>
  casadi_int CASADI_EXPORT SX::element_hash() const {
    return scalar().__hash__();
  }

  template<>
  bool CASADI_EXPORT SX::is_leaf() const {
    return scalar().is_leaf();
  }

  template<>
  bool CASADI_EXPORT SX::is_commutative() const {
    return scalar().is_commutative();
  }

  template<>
  bool CASADI_EXPORT SX::is_valid_input() const {
    for (casadi_int k=0; k<nnz(); ++k) // loop over non-zero elements
      if (!nonzeros().at(k)->is_symbolic()) // if an element is not symbolic
        return false;

    return true;
  }

  template<>
  bool CASADI_EXPORT SX::is_symbolic() const {
    if (is_dense()) {
      return is_valid_input();
    } else {
      return false;
    }
  }

  template<>
  casadi_int CASADI_EXPORT SX::op() const {
    return scalar().op();
  }

  template<>
  bool CASADI_EXPORT SX::is_op(casadi_int op) const {
    return scalar().is_op(op);
  }

  template<> bool CASADI_EXPORT SX::has_duplicates() const {
    bool has_duplicates = false;
    for (auto&& i : nonzeros_) {
      bool is_duplicate = i.get_temp()!=0;
      if (is_duplicate) {
        casadi_warning("Duplicate expression: " + str(i));
      }
      has_duplicates = has_duplicates || is_duplicate;
      i.set_temp(1);
    }
    return has_duplicates;
  }

  template<> void CASADI_EXPORT SX::reset_input() const {
    for (auto&& i : nonzeros_) {
      i.set_temp(0);
    }
  }

  template<>
  string CASADI_EXPORT SX::name() const {
    return scalar().name();
  }

  template<>
  SX CASADI_EXPORT SX::dep(casadi_int ch) const {
    return scalar().dep(ch);
  }

  template<>
  casadi_int CASADI_EXPORT SX::n_dep() const {
    return scalar().n_dep();
  }

  template<>
  void CASADI_EXPORT SX::expand(const SX& ex2, SX& ww, SX& tt) {
    casadi_assert_dev(ex2.is_scalar());
    SXElem ex = ex2.scalar();

    // Terms, weights and indices of the nodes that are already expanded
    vector<vector<SXNode*> > terms;
    vector<vector<double> > weights;
    std::map<SXNode*, casadi_int> indices;

    // Stack of nodes that are not yet expanded
    stack<SXNode*> to_be_expanded;
    to_be_expanded.push(ex.get());

    while (!to_be_expanded.empty()) { // as long as there are nodes to be expanded

      // Check if the last element on the stack is already expanded
      if (indices.find(to_be_expanded.top()) != indices.end()) {
        // Remove from stack
        to_be_expanded.pop();
        continue;
      }

      // Weights and terms
      vector<double> w; // weights
      vector<SXNode*> f; // terms

      if (to_be_expanded.top()->is_constant()) { // constant nodes are seen as multiples of one
        w.push_back(to_be_expanded.top()->to_double());
        f.push_back(casadi_limits<SXElem>::one.get());
      } else if (to_be_expanded.top()->is_symbolic()) {
        // symbolic nodes have weight one and itself as factor
        w.push_back(1);
        f.push_back(to_be_expanded.top());
      } else { // unary or binary node

        casadi_assert_dev(to_be_expanded.top()->n_dep()); // make sure that the node is binary

        // Check if addition, subtracton or multiplication
        SXNode* node = to_be_expanded.top();
        // If we have a binary node that we can factorize
        if (node->op() == OP_ADD || node->op() == OP_SUB ||
           (node->op() == OP_MUL  && (node->dep(0)->is_constant() ||
                                         node->dep(1)->is_constant()))) {
          // Make sure that both children are factorized, if not - add to stack
          if (indices.find(node->dep(0).get()) == indices.end()) {
            to_be_expanded.push(node->dep(0).get());
            continue;
          }
          if (indices.find(node->dep(1).get()) == indices.end()) {
            to_be_expanded.push(node->dep(1).get());
            continue;
          }

          // Get indices of children
          casadi_int ind1 = indices[node->dep(0).get()];
          casadi_int ind2 = indices[node->dep(1).get()];

          // If multiplication
          if (node->op() == OP_MUL) {
            double fac;
            // Multiplication where the first factor is a constant
            if (node->dep(0)->is_constant()) {
              fac = node->dep(0)->to_double();
              f = terms[ind2];
              w = weights[ind2];
            } else { // Multiplication where the second factor is a constant
              fac = node->dep(1)->to_double();
              f = terms[ind1];
              w = weights[ind1];
            }
            for (casadi_int i=0; i<w.size(); ++i) w[i] *= fac;

          } else { // if addition or subtraction
            if (node->op() == OP_ADD) {          // Addition: join both sums
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];    w.insert(w.end(), weights[ind2].begin(), weights[ind2].end());
            } else {      // Subtraction: join both sums with negative weights for second term
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];
              w.reserve(f.size());
              for (casadi_int i=0; i<weights[ind2].size(); ++i) w.push_back(-weights[ind2][i]);
            }
            // Eliminate multiple elements
            vector<double> w_new; w_new.reserve(w.size());   // weights
            vector<SXNode*> f_new;  f_new.reserve(f.size());   // terms
            std::map<SXNode*, casadi_int> f_ind; // index in f_new

            for (casadi_int i=0; i<w.size(); i++) {
              // Try to locate the node
              auto it = f_ind.find(f[i]);
              if (it == f_ind.end()) { // if the term wasn't found
                w_new.push_back(w[i]);
                f_new.push_back(f[i]);
                f_ind[f[i]] = f_new.size()-1;
              } else { // if the term already exists
                w_new[it->second] += w[i]; // just add the weight
              }
            }
            w = w_new;
            f = f_new;
          }
        } else { // if we have a binary node that we cannot factorize
          // By default,
          w.push_back(1);
          f.push_back(node);

        }
      }

      // Save factorization of the node
      weights.push_back(w);
      terms.push_back(f);
      indices[to_be_expanded.top()] = terms.size()-1;

      // Remove node from stack
      to_be_expanded.pop();
    }

    // Save expansion to output
    casadi_int thisind = indices[ex.get()];
    ww = SX(weights[thisind]);

    vector<SXElem> termsv(terms[thisind].size());
    for (casadi_int i=0; i<termsv.size(); ++i)
      termsv[i] = SXElem::create(terms[thisind][i]);
    tt = SX(termsv);
  }

  template<>
  SX CASADI_EXPORT SX::pw_const(const SX& t, const SX& tval, const SX& val) {
    // number of intervals
    casadi_int n = val.numel();

    casadi_assert(t.is_scalar(), "t must be a scalar");
    casadi_assert(tval.numel() == n-1, "dimensions do not match");

    SX ret = val->at(0);
    for (casadi_int i=0; i<n-1; ++i) {
      ret += (val(i+1)-val(i)) * (t>=tval(i));
    }

    return ret;
  }

  template<>
  SX CASADI_EXPORT SX::pw_lin(const SX& t, const SX& tval, const SX& val) {
    // Number of points
    casadi_int N = tval.numel();
    casadi_assert(N>=2, "pw_lin: N>=2");
    casadi_assert(val.numel() == N, "dimensions do not match");

    // Gradient for each line segment
    SX g = SX(1, N-1);
    for (casadi_int i=0; i<N-1; ++i)
      g(i) = (val(i+1)- val(i))/(tval(i+1)-tval(i));

    // Line segments
    SX lseg = SX(1, N-1);
    for (casadi_int i=0; i<N-1; ++i)
      lseg(i) = val(i) + g(i)*(t-tval(i));

    // Return piecewise linear function
    return pw_const(t, tval(range(1, N-1)), lseg);
  }

  template<>
  SX CASADI_EXPORT SX::gauss_quadrature(const SX& f, const SX& x, const SX& a,
      const SX& b, casadi_int order, const SX& w) {
    casadi_assert(order == 5, "gauss_quadrature: order must be 5");
    casadi_assert(w.is_empty(), "gauss_quadrature: empty weights");

    // Change variables to [-1, 1]
    if (!is_equal(a.scalar(), -1) || !is_equal(b.scalar(), 1)) {
      SX q1 = (b-a)/2;
      SX q2 = (b+a)/2;

      Function fcn("gauss_quadrature", {x}, {f});

      return q1*gauss_quadrature(fcn(q1*x+q2).at(0), x, -1, 1);
    }

    // Gauss points
    vector<double> xi;
    xi.push_back(-std::sqrt(5 + 2*std::sqrt(10.0/7))/3);
    xi.push_back(-std::sqrt(5 - 2*std::sqrt(10.0/7))/3);
    xi.push_back(0);
    xi.push_back(std::sqrt(5 - 2*std::sqrt(10.0/7))/3);
    xi.push_back(std::sqrt(5 + 2*std::sqrt(10.0/7))/3);

    // Gauss weights
    vector<double> wi;
    wi.push_back((322-13*std::sqrt(70.0))/900.0);
    wi.push_back((322+13*std::sqrt(70.0))/900.0);
    wi.push_back(128/225.0);
    wi.push_back((322+13*std::sqrt(70.0))/900.0);
    wi.push_back((322-13*std::sqrt(70.0))/900.0);

    // Evaluate at the Gauss points
    Function fcn("gauss_quadrature", {x}, {f});
    vector<SXElem> f_val(5);
    for (casadi_int i=0; i<5; ++i)
      f_val[i] = fcn(SX(xi[i])).at(0).scalar();

    // Weighted sum
    SXElem sum;
    for (casadi_int i=0; i<5; ++i)
      sum += wi[i]*f_val[i];

    return sum;
  }

  template<>
  SX CASADI_EXPORT SX::simplify(const SX& x) {
    SX r = x;
    for (casadi_int el=0; el<r.nnz(); ++el) {
      // Start by expanding the node to a weighted sum
      SX terms, weights;
      expand(r.nz(el), weights, terms);

      // Make a scalar product to get the simplified expression
      r.nz(el) = mtimes(terms.T(), weights);
    }
    return r;
  }

  template<>
  vector<SX> CASADI_EXPORT
  SX::substitute(const vector<SX>& ex, const vector<SX>& v, const vector<SX>& vdef) {

    // Assert consistent dimensions
    if (v.size()!=vdef.size()) {
      casadi_warning("subtitute: number of symbols to replace ( " + str(v.size()) + ") "
                     "must match number of expressions (" + str(vdef.size()) + ") "
                     "to replace them with.");
    }

    // Quick return if all equal
    bool all_equal = true;
    for (casadi_int k=0; k<v.size(); ++k) {
      if (v[k].size()!=vdef[k].size() || !is_equal(v[k], vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Check sparsities
    for (casadi_int k=0; k<v.size(); ++k) {
      if (v[k].sparsity()!=vdef[k].sparsity()) {
        // Expand vdef to sparsity of v if vdef is scalar
        if (vdef[k].is_scalar() && vdef[k].nnz()==1) {
          vector<SX> vdef_mod = vdef;
          vdef_mod[k] = SX(v[k].sparsity(), vdef[k]->at(0), false);
          return substitute(ex, v, vdef_mod);
        } else {
          casadi_error("Sparsities of v and vdef must match. Got v: "
                       + v[k].dim() + " and vdef: " + vdef[k].dim() + ".");
        }
      }
    }


    // Otherwise, evaluate symbolically
    Function F("tmp", v, ex);
    return F(vdef);
  }

  template<>
  SX CASADI_EXPORT SX::substitute(const SX& ex, const SX& v, const SX& vdef) {
    return substitute(vector<SX>{ex}, vector<SX>{v}, vector<SX>{vdef}).front();
  }

  template<>
  void CASADI_EXPORT SX::substitute_inplace(const vector<SX >& v, vector<SX >& vdef,
                             vector<SX >& ex, bool reverse) {
    // Assert correctness
    casadi_assert_dev(v.size()==vdef.size());
    for (casadi_int i=0; i<v.size(); ++i) {
      casadi_assert(v[i].is_symbolic(), "the variable is not symbolic");
      casadi_assert(v[i].sparsity() == vdef[i].sparsity(), "the sparsity patterns of the "
                            "expression and its defining bexpression do not match");
    }

    // Quick return if empty or single expression
    if (v.empty()) return;

    // Function inputs
    vector<SX> f_in;
    if (!reverse) f_in.insert(f_in.end(), v.begin(), v.end());

    // Function outputs
    vector<SX> f_out = vdef;
    f_out.insert(f_out.end(), ex.begin(), ex.end());

    // Write the mapping function
    Function f("tmp", f_in, f_out);

    // Get references to the internal data structures
    SXFunction *ff = f.get<SXFunction>();
    const vector<ScalarAtomic>& algorithm = ff->algorithm_;
    vector<SXElem> work(f.sz_w());

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=ff->operations_.begin();

    // Iterator to stack of constants
    vector<SXElem>::const_iterator c_it = ff->constants_.begin();

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = ff->free_vars_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_INPUT:
        // reverse is false, substitute out
        work[it->i0] = vdef.at(it->i1)->at(it->i2);
        break;
      case OP_OUTPUT:
        if (it->i0 < v.size()) {
          vdef.at(it->i0)->at(it->i2) = work[it->i1];
          if (reverse) {
            // Use the new variable henceforth, substitute in
            work[it->i1] = v.at(it->i0)->at(it->i2);
          }
        } else {
          // Auxillary output
          ex.at(it->i0 - v.size())->at(it->i2) = work[it->i1];
        }
        break;
      case OP_CONST:      work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work[it->i0] = *p_it++; break;
      default:
        {
          switch (it->op) {
            CASADI_MATH_FUN_BUILTIN(work[it->i1], work[it->i2], work[it->i0])
              }

          // Avoid creating duplicates
          const casadi_int depth = 2; // NOTE: a higher depth could possibly give more savings
          work[it->i0].assignIfDuplicate(*b_it++, depth);
        }
      }
    }
  }

  template<>
  bool CASADI_EXPORT SX::depends_on(const SX &x, const SX &arg) {
    if (x.nnz()==0) return false;

    // Construct a temporary algorithm
    Function temp("temp", {arg}, {x});

    // Perform a single dependency sweep
    vector<bvec_t> t_in(arg.nnz(), 1), t_out(x.nnz());
    temp({get_ptr(t_in)}, {get_ptr(t_out)});

    // Loop over results
    for (casadi_int i=0; i<t_out.size(); ++i) {
      if (t_out[i]) return true;
    }

    return false;
  }

  template<>
  SX CASADI_EXPORT SX::jacobian(const SX &f, const SX &x, const Dict& opts) {
    // Propagate verbose option to helper function
    Dict h_opts;
    Dict opts_remainder = extract_from_dict(opts, "helper_options", h_opts);
    Function h("jac_helper", {x}, {f}, h_opts);
    return h.get<SXFunction>()->jac(0, 0, opts_remainder);
  }

  template<>
  SX CASADI_EXPORT SX::hessian(const SX &ex, const SX &arg, SX &g, const Dict& opts) {
    Dict all_opts = opts;
    if (!opts.count("symmetric")) all_opts["symmetric"] = true;
    g = gradient(ex, arg);
    return jacobian(g, arg, all_opts);
  }

  template<>
  SX CASADI_EXPORT SX::hessian(const SX &ex, const SX &arg, const Dict& opts) {
    SX g;
    return hessian(ex, arg, g, opts);
  }

  template<>
  std::vector<std::vector<SX> > CASADI_EXPORT
  SX::forward(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts) {

    Dict h_opts;
    Dict opts_remainder = extract_from_dict(opts, "helper_options", h_opts);
    // Read options
    bool always_inline = false;
    bool never_inline = false;
    for (auto&& op : opts_remainder) {
      if (op.first=="always_inline") {
        always_inline = op.second;
      } else if (op.first=="never_inline") {
        never_inline = op.second;
      } else {
        casadi_error("No such option: " + string(op.first));
      }
    }
    // Call internal function on a temporary object
    Function temp("forward_temp", arg, ex, h_opts);
    std::vector<std::vector<SX> > ret;
    temp->call_forward(arg, ex, v, ret, always_inline, never_inline);
    return ret;
  }

  template<>
  std::vector<std::vector<SX> > CASADI_EXPORT
  SX::reverse(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts) {

    Dict h_opts;
    Dict opts_remainder = extract_from_dict(opts, "helper_options", h_opts);

    // Read options
    bool always_inline = false;
    bool never_inline = false;
    for (auto&& op : opts_remainder) {
      if (op.first=="always_inline") {
        always_inline = op.second;
      } else if (op.first=="never_inline") {
        never_inline = op.second;
      } else {
        casadi_error("No such option: " + string(op.first));
      }
    }
    // Call internal function on a temporary object
    Function temp("reverse_temp", arg, ex, h_opts);
    std::vector<std::vector<SX> > ret;
    temp->call_reverse(arg, ex, v, ret, always_inline, never_inline);
    return ret;
  }

  template<>
  std::vector<bool> CASADI_EXPORT SX::which_depends(const SX &expr,
      const SX &var, casadi_int order, bool tr) {
    return _which_depends(expr, var, order, tr);
  }

  template<>
  SX CASADI_EXPORT SX::taylor(const SX& f, const SX& x,
                const SX& a, casadi_int order) {
    casadi_assert_dev(x.is_scalar() && a.is_scalar());
    if (f.nnz()!=f.numel())
      throw CasadiException("taylor: not implemented for sparse matrices");
    SX ff = vec(f.T());

    SX result = substitute(ff, x, a);
    double nf=1;
    SX dx = (x-a);
    SX dxa = (x-a);
    for (casadi_int i=1; i<=order; i++) {
      ff = jacobian(ff, x);
      nf*=static_cast<double>(i);
      result+=1/nf * substitute(ff, x, a) * dxa;
      dxa*=dx;
    }
    return reshape(result, f.size2(), f.size1()).T();
  }

  SX mtaylor_recursive(const SX& ex, const SX& x, const SX& a, casadi_int order,
                       const vector<casadi_int>&order_contributions,
                       const SXElem & current_dx=casadi_limits<SXElem>::one,
                       double current_denom=1, casadi_int current_order=1) {
    SX result = substitute(ex, x, a)*current_dx/current_denom;
    for (casadi_int i=0;i<x.nnz();i++) {
      if (order_contributions[i]<=order) {
        result += mtaylor_recursive(SX::jacobian(ex, x->at(i)),
                                    x, a,
                                    order-order_contributions[i],
                                    order_contributions,
                                    current_dx*(x->at(i)-a->at(i)),
                                    current_denom*static_cast<double>(current_order),
                                    current_order+1);
      }
    }
    return result;
  }

  template<>
  SX CASADI_EXPORT SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order,
                 const vector<casadi_int>& order_contributions) {
    casadi_assert(f.nnz()==f.numel() && x.nnz()==x.numel(),
                          "mtaylor: not implemented for sparse matrices");

    casadi_assert(x.nnz()==order_contributions.size(),
                          "mtaylor: number of non-zero elements in x (" + str(x.nnz())
                          + ") must match size of order_contributions ("
                          + str(order_contributions.size()) + ")");

    return reshape(mtaylor_recursive(vec(f), x, a, order,
                                     order_contributions),
                   f.size2(), f.size1()).T();
  }

  template<>
  SX CASADI_EXPORT SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order) {
    return mtaylor(f, x, a, order, vector<casadi_int>(x.nnz(), 1));
  }

  template<>
  casadi_int CASADI_EXPORT SX::n_nodes(const SX& x) {
    Function f("tmp", {SX()}, {x});
    return f.n_nodes();
  }

  template<>
  string CASADI_EXPORT
  SX::print_operator(const SX& X, const vector<string>& args) {
    SXElem x = X.scalar();
    casadi_int ndeps = casadi_math<double>::ndeps(x.op());
    casadi_assert(ndeps==1 || ndeps==2, "Not a unary or binary operator");
    casadi_assert(args.size()==ndeps, "Wrong number of arguments");
    if (ndeps==1) {
      return casadi_math<double>::print(x.op(), args.at(0));
    } else {
      return casadi_math<double>::print(x.op(), args.at(0), args.at(1));
    }
  }

  template<>
  vector<SX> CASADI_EXPORT SX::symvar(const SX& x) {
    Function f("tmp", vector<SX>{}, {x});
    return f.free_sx();
  }

  template<>
  void CASADI_EXPORT SX::shared(vector<SX >& ex,
                         vector<SX >& v_sx,
                         vector<SX >& vdef_sx,
                         const string& v_prefix,
                         const string& v_suffix) {

    // Sort the expression
    Function f("tmp", vector<SX>(), ex);
    SXFunction *ff = f.get<SXFunction>();

    // Get references to the internal data structures
    const vector<ScalarAtomic>& algorithm = ff->algorithm_;
    vector<SXElem> work(f.sz_w());
    vector<SXElem> work2 = work;

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=ff->operations_.begin();

    // Iterator to stack of constants
    vector<SXElem>::const_iterator c_it = ff->constants_.begin();

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = ff->free_vars_.begin();

    // Count how many times an expression has been used
    vector<casadi_int> usecount(work.size(), 0);

    // Evaluate the algorithm
    vector<SXElem> v, vdef;
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      // Increase usage counters
      switch (it->op) {
      case OP_CONST:
      case OP_PARAMETER:
        break;
        CASADI_MATH_BINARY_BUILTIN // Binary operation
          if (usecount[it->i2]==0) {
            usecount[it->i2]=1;
          } else if (usecount[it->i2]==1) {
            // Get a suitable name
            vdef.push_back(work[it->i2]);
            usecount[it->i2]=-1; // Extracted, do not extract again
          }
        // fall-through
      case OP_OUTPUT:
      default: // Unary operation, binary operation or output
        if (usecount[it->i1]==0) {
          usecount[it->i1]=1;
        } else if (usecount[it->i1]==1) {
          vdef.push_back(work[it->i1]);
          usecount[it->i1]=-1; // Extracted, do not extract again
        }
      }

      // Perform the operation
      switch (it->op) {
      case OP_OUTPUT:
        break;
      case OP_CONST:
      case OP_PARAMETER:
        usecount[it->i0] = -1; // Never extract since it is a primitive type
        break;
      default:
        work[it->i0] = *b_it++;
        usecount[it->i0] = 0; // Not (yet) extracted
        break;
      }
    }

    // Create intermediate variables
    stringstream v_name;
    for (casadi_int i=0; i<vdef.size(); ++i) {
      v_name.str(string());
      v_name << v_prefix << i << v_suffix;
      v.push_back(SXElem::sym(v_name.str()));
    }

    casadi_assert(vdef.size()<numeric_limits<int>::max(), "Integer overflow");

    // Mark the above expressions
    for (casadi_int i=0; i<vdef.size(); ++i) {
      vdef[i].set_temp(static_cast<int>(i)+1);
    }

    // Save the marked nodes for later cleanup
    vector<SXElem> marked = vdef;

    // Reset iterator
    b_it=ff->operations_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_OUTPUT:     ex.at(it->i0)->at(it->i2) = work[it->i1];      break;
      case OP_CONST:      work2[it->i0] = work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work2[it->i0] = work[it->i0] = *p_it++; break;
      default:
        {
          switch (it->op) {
            CASADI_MATH_FUN_BUILTIN(work[it->i1], work[it->i2], work[it->i0])
              }
          work2[it->i0] = *b_it++;

          // Replace with intermediate variables
          casadi_int ind = work2[it->i0].get_temp()-1;
          if (ind>=0) {
            vdef.at(ind) = work[it->i0];
            work[it->i0] = v.at(ind);
          }
        }
      }
    }

    // Unmark the expressions
    for (vector<SXElem>::iterator it=marked.begin(); it!=marked.end(); ++it) {
      it->set_temp(0);
    }

    // Save v, vdef
    v_sx.resize(v.size());
    copy(v.begin(), v.end(), v_sx.begin());
    vdef_sx.resize(vdef.size());
    copy(vdef.begin(), vdef.end(), vdef_sx.begin());
  }

  template<>
  SX CASADI_EXPORT SX::poly_coeff(const SX& ex, const SX& x) {
    casadi_assert_dev(ex.is_scalar());
    casadi_assert_dev(x.is_scalar());
    casadi_assert_dev(x.is_symbolic());

    vector<SXElem> r;

    SX j = ex;
    casadi_int mult = 1;
    bool success = false;
    for (casadi_int i=0; i<1000; ++i) {
      r.push_back((substitute(j, x, 0)/static_cast<double>(mult)).scalar());
      j = jacobian(j, x);
      if (j.nnz()==0) {
        success = true;
        break;
      }
      mult*=i+1;
    }

    if (!success) casadi_error("poly: supplied expression does not appear to be polynomial.");

    std::reverse(r.begin(), r.end());

    return r;
  }

  template<>
  SX CASADI_EXPORT SX::poly_roots(const SX& p) {
    casadi_assert(p.size2()==1,
                          "poly_root(): supplied parameter must be column vector but got "
                          + p.dim() + ".");
    casadi_assert_dev(p.is_dense());
    if (p.size1()==2) { // a*x + b
      SX a = p(0);
      SX b = p(1);
      return -b/a;
    } else if (p.size1()==3) { // a*x^2 + b*x + c
      SX a = p(0);
      SX b = p(1);
      SX c = p(2);
      SX ds = sqrt(b*b-4*a*c);
      SX bm = -b;
      SX a2 = 2*a;
      SX ret = SX::vertcat({(bm-ds)/a2, (bm+ds)/a2});
      return ret;
    } else if (p.size1()==4) {
      // www.cs.iastate.edu/~cs577/handouts/polyroots.pdf
      SX ai = 1/p(0);

      SX p_ = p(1)*ai;
      SX q  = p(2)*ai;
      SX r  = p(3)*ai;

      SX pp = p_*p_;

      SX a = q - pp/3;
      SX b = r + 2.0/27*pp*p_-p_*q/3;

      SX a3 = a/3;

      SX phi = acos(-b/2/sqrt(-a3*a3*a3));

      SX ret = SX::vertcat({cos(phi/3), cos((phi+2*pi)/3), cos((phi+4*pi)/3)});
      ret*= 2*sqrt(-a3);

      ret-= p_/3;
      return ret;
    } else if (p.size1()==5) {
      SX ai = 1/p(0);
      SX b = p(1)*ai;
      SX c = p(2)*ai;
      SX d = p(3)*ai;
      SX e = p(4)*ai;

      SX bb= b*b;
      SX f = c - (3*bb/8);
      SX g = d + (bb*b / 8) - b*c/2;
      SX h = e - (3*bb*bb/256) + (bb * c/16) - (b*d/4);
      SX poly = SX::vertcat({1, f/2, ((f*f -4*h)/16), -g*g/64});
      SX y = poly_roots(poly);

      SX r0 = y(0); // NOLINT(cppcoreguidelines-slicing)
      SX r1 = y(2); // NOLINT(cppcoreguidelines-slicing)

      SX p = sqrt(r0); // two non-zero-roots
      SX q = sqrt(r1);

      SX r = -g/(8*p*q);

      SX s = b/4;

      SX ret = SX::vertcat({
          p + q + r -s,
            p - q - r -s,
            -p + q - r -s,
            -p - q + r -s});
      return ret;
    } else if (is_equal(p(p.nnz()-1)->at(0), 0)) {
      SX ret = SX::vertcat({poly_roots(p(range(p.nnz()-1))), 0});
      return ret;
    } else {
      casadi_error("poly_root(): can only solve cases for first or second order polynomial. "
                   "Got order " + str(p.size1()-1) + ".");
    }

  }

  template<>
  SX CASADI_EXPORT SX::eig_symbolic(const SX& m) {
    casadi_assert(m.size1()==m.size2(), "eig(): supplied matrix must be square");

    vector<SX> ret;

    /// Bring m in block diagonal form, calculating eigenvalues of each block separately
    vector<casadi_int> offset;
    vector<casadi_int> index;
    casadi_int nb = m.sparsity().scc(offset, index);

    SX m_perm = m(offset, offset);

    SX l = SX::sym("l");

    for (casadi_int k=0; k<nb; ++k) {
      vector<casadi_int> r = range(index.at(k), index.at(k+1));
      // det(lambda*I-m) = 0
      ret.push_back(poly_roots(poly_coeff(det(SX::eye(r.size())*l-m_perm(r, r)), l)));
    }

    return vertcat(ret);
  }

  template<>
  void CASADI_EXPORT SX::print_split(casadi_int nnz, const SXElem* nonzeros, vector<string>& nz,
                      vector<string>& inter) {
    // Find out which noded can be inlined
    std::map<const SXNode*, casadi_int> nodeind;
    for (casadi_int i=0; i<nnz; ++i) nonzeros[i]->can_inline(nodeind);

    // Print expression
    nz.resize(0);
    nz.reserve(nnz);
    inter.resize(0);
    for (casadi_int i=0; i<nnz; ++i) nz.push_back(nonzeros[i]->print_compact(nodeind, inter));
  }

  template<> vector<SX> CASADI_EXPORT SX::get_input(const Function& f) {
    return f.sx_in();
  }

  template<> vector<SX> CASADI_EXPORT SX::get_free(const Function& f) {
    return f.free_sx();
  }

  template<>
  Dict CASADI_EXPORT SX::info() const {
    return {{"function", Function("f", std::vector<SX>{}, std::vector<SX>{*this})}};
  }

  template<>
  void CASADI_EXPORT SX::to_file(const std::string& filename,
      const Sparsity& sp, const SXElem* nonzeros,
      const std::string& format_hint) {
    casadi_error("Not implemented");
  }

  template class CASADI_EXPORT Matrix< SXElem >;
} // namespace casadi
