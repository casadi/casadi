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


#include "matrix_impl.hpp"
#include "../function/linear_solver.hpp"
#include "../function/sx_function_internal.hpp"
#include "../sx/sx_node.hpp"

using namespace std;

namespace casadi {

  template<>
  bool Matrix<int>::isSlice(bool ind1) const {
    return isscalar() || (iscolumn() && isdense() && Slice::isSlice(data(), ind1));
  }

  template<>
  Slice Matrix<int>::toSlice(bool ind1) const {
    return isscalar() ? Slice(toScalar(), ind1) : Slice(data(), ind1);
  }

  template<>
  Matrix<double> Matrix<double>::
  zz_solve(const Matrix<double>& b,
           const std::string& lsolver, const Dict& dict) const {
    const Matrix<double>& A = *this;
    LinearSolver mysolver("tmp", lsolver, A.sparsity(), b.size2(), dict);
    mysolver.setInput(A, LINSOL_A);
    mysolver.setInput(b, LINSOL_B);
    mysolver.prepare();
    mysolver.solve(false);

    return mysolver.output(LINSOL_X);
  }

  template<>
  Matrix<double> Matrix<double>::
  zz_pinv(const std::string& lsolver,
          const Dict& dict) const {
    const Matrix<double>& A = *this;
    if (A.size1()>=A.size2()) {
      return solve(mul(A.T(), A), A.T(), lsolver, dict);
    } else {
      return solve(mul(A, A.T()), A, lsolver, dict).T();
    }
  }

  template<>
  bool Matrix<SXElement>::__nonzero__() const {
    if (numel()!=1) {casadi_error("Only scalar Matrix could have a truth value, but you "
                                  "provided a shape" << dimString());}
    return at(0).__nonzero__();
  }

  template<>
  void SX::setEqualityCheckingDepth(int eq_depth) {
    SXNode::eq_depth_ = eq_depth;
  }

  template<>
  int SX::getEqualityCheckingDepth() {
    return SXNode::eq_depth_;
  }

  template<>
  SX GenericMatrix<SX>::sym(const std::string& name, const Sparsity& sp) {
    // Create a dense n-by-m matrix
    std::vector<SXElement> retv;

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
          retv.push_back(SXElement::sym(varname));
      }
    } else if (sp.isscalar(true)) {
      retv.push_back(SXElement::sym(name));
    } else {
      // Scalar
      std::stringstream ss;
      for (int k=0; k<sp.nnz(); ++k) {
        ss.str("");
        ss << name << "_" << k;
        retv.push_back(SXElement::sym(ss.str()));
      }
    }

    // Determine dimensions automatically if empty
    if (sp.isscalar(true)) {
      return SX(retv);
    } else {
      return SX(sp, retv, false);
    }
  }

  template<>
  bool SX::isRegular() const {
    // First pass: ignore symbolics
    for (int i=0; i<nnz(); ++i) {
      const SXElement& x = at(i);
      if (x.isConstant()) {
        if (x.isNan() || x.isInf() || x.isMinusInf()) return false;
      }
    }
    // Second pass: don't ignore symbolics
    for (int i=0; i<nnz(); ++i) {
      if (!at(i).isRegular()) return false;
    }
    return true;
  }

  template<>
  bool SX::isSmooth() const {
    // Make a function
    SXFunction temp("temp", make_vector(SX()), make_vector(*this));

    // Run the function on the temporary variable
    return temp->isSmooth();
  }

  template<>
  size_t SX::getElementHash() const {
    return toScalar().__hash__();
  }

  template<>
  bool SX::isLeaf() const {
    return toScalar().isLeaf();
  }

  template<>
  bool SX::isCommutative() const {
    return toScalar().isCommutative();
  }

  template<>
  bool SX::isSymbolic() const {
    if (isdense()) {
      return isValidInput();
    } else {
      return false;
    }
  }

  template<>
  bool SX::isValidInput() const {
    for (int k=0; k<nnz(); ++k) // loop over non-zero elements
      if (!at(k)->isSymbolic()) // if an element is not symbolic
        return false;

    return true;
  }

  template<> bool SX::hasDuplicates() {
    bool has_duplicates = false;
    for (vector<SXElement>::iterator it = begin(); it != end(); ++it) {
      bool is_duplicate = it->getTemp()!=0;
      if (is_duplicate) {
        userOut<true, PL_WARN>() << "Duplicate expression: " << *it << endl;
      }
      has_duplicates = has_duplicates || is_duplicate;
      it->setTemp(1);
    }
    return has_duplicates;
  }

  template<> void SX::resetInput() {
    for (vector<SXElement>::iterator it = begin(); it != end(); ++it) {
      it->setTemp(0);
    }
  }

  template<>
  double SX::getValue(int k) const {
    return at(k).getValue();
  }

  template<>
  int SX::getIntValue() const {
    return toScalar().getIntValue();
  }

  template<>
  std::vector<double> SX::nonzeros() const {
    std::vector<double> ret(nnz());
    for (size_t i=0; i<ret.size(); ++i) {
      ret[i] = at(i).getValue();
    }
    return ret;
  }

  template<>
  std::vector<int> SX::nonzeros_int() const {
    std::vector<int> ret(nnz());
    for (size_t i=0; i<ret.size(); ++i) {
      ret[i] = at(i).getIntValue();
    }
    return ret;
  }

  template<>
  std::string SX::getName() const {
    return toScalar().getName();
  }

  template<>
  SX SX::getDep(int ch) const {
    return toScalar().getDep(ch);
  }

  template<>
  int SX::getNdeps() const {
    return toScalar().getNdeps();
  }

  template<>
  void SX::zz_expand(SX &ww, SX& tt) const {
    const SX& ex2 = *this;
    casadi_assert(ex2.isscalar());
    SXElement ex = ex2.toScalar();

    // Terms, weights and indices of the nodes that are already expanded
    std::vector<std::vector<SXNode*> > terms;
    std::vector<std::vector<double> > weights;
    std::map<SXNode*, int> indices;

    // Stack of nodes that are not yet expanded
    std::stack<SXNode*> to_be_expanded;
    to_be_expanded.push(ex.get());

    while (!to_be_expanded.empty()) { // as long as there are nodes to be expanded

      // Check if the last element on the stack is already expanded
      if (indices.find(to_be_expanded.top()) != indices.end()) {
        // Remove from stack
        to_be_expanded.pop();
        continue;
      }

      // Weights and terms
      std::vector<double> w; // weights
      std::vector<SXNode*> f; // terms

      if (to_be_expanded.top()->isConstant()) { // constant nodes are seen as multiples of one
        w.push_back(to_be_expanded.top()->getValue());
        f.push_back(casadi_limits<SXElement>::one.get());
      } else if (to_be_expanded.top()->isSymbolic()) {
        // symbolic nodes have weight one and itself as factor
        w.push_back(1);
        f.push_back(to_be_expanded.top());
      } else { // binary node

        casadi_assert(to_be_expanded.top()->hasDep()); // make sure that the node is binary

        // Check if addition, subtracton or multiplication
        SXNode* node = to_be_expanded.top();
        // If we have a binary node that we can factorize
        if (node->getOp() == OP_ADD || node->getOp() == OP_SUB ||
           (node->getOp() == OP_MUL  && (node->dep(0)->isConstant() ||
                                         node->dep(1)->isConstant()))) {
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
          int ind1 = indices[node->dep(0).get()];
          int ind2 = indices[node->dep(1).get()];

          // If multiplication
          if (node->getOp() == OP_MUL) {
            double fac;
            if (node->dep(0)->isConstant()) { // Multiplication where the first factor is a constant
              fac = node->dep(0)->getValue();
              f = terms[ind2];
              w = weights[ind2];
            } else { // Multiplication where the second factor is a constant
              fac = node->dep(1)->getValue();
              f = terms[ind1];
              w = weights[ind1];
            }
            for (int i=0; i<w.size(); ++i) w[i] *= fac;

          } else { // if addition or subtraction
            if (node->getOp() == OP_ADD) {          // Addition: join both sums
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];    w.insert(w.end(), weights[ind2].begin(), weights[ind2].end());
            } else {      // Subtraction: join both sums with negative weights for second term
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];
              w.reserve(f.size());
              for (int i=0; i<weights[ind2].size(); ++i) w.push_back(-weights[ind2][i]);
            }
            // Eliminate multiple elements
            std::vector<double> w_new; w_new.reserve(w.size());   // weights
            std::vector<SXNode*> f_new;  f_new.reserve(f.size());   // terms
            std::map<SXNode*, int> f_ind; // index in f_new

            for (int i=0; i<w.size(); i++) {
              // Try to locate the node
              std::map<SXNode*, int>::iterator it = f_ind.find(f[i]);
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
    int thisind = indices[ex.get()];
    ww = SX(weights[thisind]);

    vector<SXElement> termsv(terms[thisind].size());
    for (int i=0; i<termsv.size(); ++i)
      termsv[i] = SXElement::create(terms[thisind][i]);
    tt = SX(termsv);
  }

  template<>
  SX SX::zz_pw_const(const SX &tval,
                     const SX &val) const {
    const SX &t = *this;
    // number of intervals
    int n = val.numel();

    casadi_assert_message(t.isscalar(), "t must be a scalar");
    casadi_assert_message(tval.numel() == n-1, "dimensions do not match");

    SX ret = val.at(0);
    for (int i=0; i<n-1; ++i) {
      ret += (val(0, i+1)-val(0, i)) * (t>=tval(0, i));
    }

    return ret;
  }

  template<>
  SX SX::zz_pw_lin(const SX &tval,
                   const SX &val) const {
    const SX &t = *this;

    // Number of points
    int N = tval.numel();
    casadi_assert_message(N>=2, "pw_lin: N>=2");

    // Gradient for each line segment
    SX g = SX(1, N-1);
    for (int i=0; i<N-1; ++i) {
      g(0, i) = (val(0, i+1)- val(0, i))/(tval(0, i+1)-tval(0, i));
    }

    // Line segments
    SX lseg = SX(1, N-1);
    for (int i=0; i<N-1; ++i)
      lseg(0, i) = val(0, i) + g(0, i)*(t-tval(0, i));

    // interior time points
    SX tint = tval(0, range(N-2));

    // Return piecewise linear function
    return pw_const(t, tint, lseg);
  }

  template<>
  SX SX::zz_if_else(const SX &if_true, const SX &if_false, bool short_circuit) const {
    return if_else_zero(*this, if_true) + if_else_zero(!*this, if_false);
  }

  template<>
  SX SX::zz_gauss_quadrature(const SX &x,
                             const SX &a,
                             const SX &b, int order,
                             const SX& w) const {
    const SX &f = *this;
    casadi_assert_message(order == 5, "gauss_quadrature: order must be 5");
    casadi_assert_message(w.isempty(), "gauss_quadrature: empty weights");

    // Change variables to [-1, 1]
    if (!isEqual(a.toScalar(), -1) || !isEqual(b.toScalar(), 1)) {
      SX q1 = (b-a)/2;
      SX q2 = (b+a)/2;

      SXFunction fcn("gauss_quadrature", make_vector(x), make_vector(f));

      return q1*gauss_quadrature(fcn(q1*x+q2).at(0), x, -1, 1);
    }

    // Gauss points
    vector<double> xi;
    xi.push_back(-sqrt(5 + 2*sqrt(10.0/7))/3);
    xi.push_back(-sqrt(5 - 2*sqrt(10.0/7))/3);
    xi.push_back(0);
    xi.push_back(sqrt(5 - 2*sqrt(10.0/7))/3);
    xi.push_back(sqrt(5 + 2*sqrt(10.0/7))/3);

    // Gauss weights
    vector<double> wi;
    wi.push_back((322-13*sqrt(70.0))/900.0);
    wi.push_back((322+13*sqrt(70.0))/900.0);
    wi.push_back(128/225.0);
    wi.push_back((322+13*sqrt(70.0))/900.0);
    wi.push_back((322-13*sqrt(70.0))/900.0);

    // Evaluate at the Gauss points
    SXFunction fcn("gauss_quadrature", make_vector(x), make_vector(f));
    vector<SXElement> f_val(5);
    for (int i=0; i<5; ++i)
      f_val[i] = fcn(SX(xi[i])).at(0).toScalar();

    // Weighted sum
    SXElement sum;
    for (int i=0; i<5; ++i)
      sum += wi[i]*f_val[i];

    return sum;
  }

  template<>
  SX SX::zz_simplify() const {
    SX ex = *this;
    for (int el=0; el<ex.nnz(); ++el) ex.at(el) = simplify(ex.at(el));
    return ex;
  }

  template<>
  SX SX::zz_substitute(const SX& v,
                       const SX& vdef) const {
    return substitute(vector<SX>(1, *this), vector<SX>(1, v), vector<SX>(1, vdef)).front();
  }

  template<>
  std::vector<SX >
  SX::zz_substitute(const std::vector<SX >& ex,
                    const std::vector<SX >& v,
                    const std::vector<SX >& vdef) {

    // Assert consistent dimensions
    casadi_assert_warning(v.size()==vdef.size(), "subtitute: number of symbols to replace ( "
                          << v.size() << ") must match number of expressions (" << vdef.size()
                          << ") to replace them with.");

    // Quick return if all equal
    bool all_equal = true;
    for (int k=0; k<v.size(); ++k) {
      if (v[k].size()!=vdef[k].size() || !isEqual(v[k], vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Check sparsities
    for (int k=0; k<v.size(); ++k) {
      if (v[k].sparsity()!=vdef[k].sparsity()) {
        // Expand vdef to sparsity of v if vdef is scalar
        if (vdef[k].isscalar() && vdef[k].nnz()==1) {
          std::vector<SX> vdef_mod = vdef;
          vdef_mod[k] = SX(v[k].sparsity(), vdef[k].at(0), false);
          return substitute(ex, v, vdef_mod);
        } else {
          casadi_error("subsitute(ex, v, vdef): sparsities of v and vdef must match. Got v: "
                       << v[k].dimString() << " and " << "vdef: " << vdef[k].dimString() << ".");
        }
      }
    }


    // Otherwise, evaluate symbolically
    SXFunction F("tmp", v, ex);
    return F(vdef);
  }

  template<>
  void SX::zz_substituteInPlace(const std::vector<SX >& v,
                                std::vector<SX >& vdef,
                                std::vector<SX >& ex,
                                bool reverse) {
    // Assert correctness
    casadi_assert(v.size()==vdef.size());
    for (int i=0; i<v.size(); ++i) {
      casadi_assert_message(v[i].isSymbolic(), "the variable is not symbolic");
      casadi_assert_message(v[i].sparsity() == vdef[i].sparsity(), "the sparsity patterns of the "
                            "expression and its defining bexpression do not match");
    }

    // Quick return if empty or single expression
    if (v.empty()) return;

    // Function inputs
    std::vector<SX> f_in;
    if (!reverse) f_in.insert(f_in.end(), v.begin(), v.end());

    // Function outputs
    std::vector<SX> f_out = vdef;
    f_out.insert(f_out.end(), ex.begin(), ex.end());

    // Write the mapping function
    SXFunction f("tmp", f_in, f_out);

    // Get references to the internal data structures
    const vector<ScalarAtomic>& algorithm = f.algorithm();
    vector<SXElement> work(f.getWorkSize());

    // Iterator to the binary operations
    vector<SXElement>::const_iterator b_it=f->operations_.begin();

    // Iterator to stack of constants
    vector<SXElement>::const_iterator c_it = f->constants_.begin();

    // Iterator to free variables
    vector<SXElement>::const_iterator p_it = f->free_vars_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_INPUT:
        // reverse is false, substitute out
        work[it->i0] = vdef.at(it->i1).at(it->i2);
        break;
      case OP_OUTPUT:
        if (it->i0 < v.size()) {
          vdef.at(it->i0).at(it->i2) = work[it->i1];
          if (reverse) {
            // Use the new variable henceforth, substitute in
            work[it->i1] = v.at(it->i0).at(it->i2);
          }
        } else {
          // Auxillary output
          ex.at(it->i0 - v.size()).at(it->i2) = work[it->i1];
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
          const int depth = 2; // NOTE: a higher depth could possibly give more savings
          work[it->i0].assignIfDuplicate(*b_it++, depth);
        }
      }
    }
  }

  template<>
  bool SX::zz_dependsOn(const SX &arg) const {
    if (nnz()==0) return false;

    // Construct a temporary algorithm
    SXFunction temp("temp", make_vector(arg), make_vector(*this));
    temp.spInit(true);

    bvec_t* input_ =  get_bvec_t(temp.input().data());
    // Make a column with all variables active
    std::fill(input_, input_+temp.input().nnz(), bvec_t(1));
    bvec_t* output_ = get_bvec_t(temp.output().data());
    // Perform a single dependency sweep
    temp.spEvaluate(true);

    // Loop over results
    for (int i=0;i<temp.output().nnz();++i) {
      if (output_[i]) return true;
    }

    return false;
  }

  template<>
  SX SX::zz_jacobian(const SX &arg) const {
    SXFunction temp("temp", make_vector(arg), make_vector(*this)); // make a runtime
    return SX::jac(temp);
  }

  template<>
  SX SX::zz_gradient(const SX &arg) const {
    SXFunction temp("temp", make_vector(arg), make_vector(*this)); // make a runtime
    return SX::grad(temp);
  }

  template<>
  SX SX::zz_tangent(const SX &arg) const {
    SXFunction temp("temp", make_vector(arg), make_vector(*this)); // make a runtime
    return SX::tang(temp);
  }

  template<>
  SX SX::zz_hessian(const SX &arg) const {
    SX g;
    return hessian(*this, arg, g);
  }

  template<>
  SX SX::zz_hessian(const SX &arg, SX &g) const {
    g = gradient(*this, arg);
    SXFunction gfcn("gfcn", make_vector(arg), make_vector(g)); // make a runtime
    return SX::jac(gfcn, 0, 0, false, true);
  }

  template<>
  SX SX::zz_jacobianTimesVector(const SX &arg, const SX &v, bool transpose_jacobian) const {
    SXFunction f("tmp", make_vector(arg), make_vector(*this));

    // Split up v
    vector<SX> vv = horzsplit(v);

    // Make sure well-posed
    casadi_assert(vv.size() >= 1);
    casadi_assert(iscolumn());
    casadi_assert(arg.iscolumn());
    if (transpose_jacobian) {
      casadi_assert(v.size1()==size1());
    } else {
      casadi_assert(v.size1()==arg.size1());
    }

    // Number of sensitivities
    int nfsens = transpose_jacobian ? 0 : vv.size();
    int nasens = transpose_jacobian ? vv.size() : 0;

    // Assemble arguments and directional derivatives
    vector<SX> argv = f.sx_in();
    vector<SX> resv = f(argv);
    vector<vector<SX> > fseed(nfsens, argv), fsens(nfsens, resv),
        aseed(nasens, resv), asens(nasens, argv);

    for (int dir=0; dir<vv.size(); ++dir) {
      if (transpose_jacobian) {
        aseed[dir][0].set(vv[dir]);
      } else {
        fseed[dir][0].set(vv[dir]);
      }
    }

    // Evaluate with directional derivatives, output is the same as the funciton inputs
    f.callDerivative(argv, resv, fseed, fsens, aseed, asens);

    // Get the results
    for (int dir=0; dir<vv.size(); ++dir) {
      if (transpose_jacobian) {
        vv[dir] = asens[dir][0];
      } else {
        vv[dir] = fsens[dir][0];
      }
    }
    return horzcat(vv);
  }

  template<>
  SX SX::zz_taylor(const SX& x,
                   const SX& a, int order) const {
    casadi_assert(x.isscalar() && a.isscalar());
    if (nnz()!=numel())
      throw CasadiException("taylor: not implemented for sparse matrices");
    SX ff = vec(T());

    SX result = substitute(ff, x, a);
    double nf=1;
    SX dx = (x-a);
    SX dxa = (x-a);
    for (int i=1;i<=order;i++) {
      ff = jacobian(ff, x);
      nf*=i;
      result+=1/nf * substitute(ff, x, a) * dxa;
      dxa*=dx;
    }
    return reshape(result, size2(), size1()).T();
  }

  template<>
  SX SX::zz_mtaylor(const SX& x, const SX& a, int order) const {
    return mtaylor(*this, x, a, order, std::vector<int>(x.nnz(), 1));
  }

  SX mtaylor_recursive(const SX& ex, const SX& x, const SX& a, int order,
                       const std::vector<int>&order_contributions,
                       const SXElement & current_dx=casadi_limits<SXElement>::one,
                       double current_denom=1, int current_order=1) {
    SX result = substitute(ex, x, a)*current_dx/current_denom;
    for (int i=0;i<x.nnz();i++) {
      if (order_contributions[i]<=order) {
        result += mtaylor_recursive(
                                    jacobian(ex, x.at(i)),
                                    x, a,
                                    order-order_contributions[i],
                                    order_contributions,
                                    current_dx*(x.at(i)-a.at(i)),
                                    current_denom*current_order, current_order+1);
      }
    }
    return result;
  }

  template<>
  SX SX::zz_mtaylor(const SX& x, const SX& a, int order,
                    const std::vector<int>& order_contributions) const {
    casadi_assert_message(nnz()==numel() && x.nnz()==x.numel(),
                          "mtaylor: not implemented for sparse matrices");

    casadi_assert_message(x.nnz()==order_contributions.size(),
                          "mtaylor: number of non-zero elements in x (" <<  x.nnz()
                          << ") must match size of order_contributions ("
                          << order_contributions.size() << ")");

    return reshape(mtaylor_recursive(vec(*this), x, a, order,
                                     order_contributions),
                   size2(), size1()).T();
  }

  template<>
  int SX::zz_countNodes() const {
    SXFunction f("tmp", make_vector(SX()), make_vector(*this));
    return f.countNodes();
  }

  template<>
  std::string
  SX::zz_getOperatorRepresentation(const std::vector<std::string>& args) const {
    SXElement x = toScalar();
    if (!x.hasDep())
        throw CasadiException("getOperatorRepresentation: SXElement must be binary operator");
    if (args.size() == 0 || (casadi_math<double>::ndeps(x.getOp())==2 && args.size() < 2))
        throw CasadiException("getOperatorRepresentation: not enough arguments supplied");
    std::stringstream s;
    casadi_math<double>::print(x.getOp(), s, args[0], args[1]);
    return s.str();
  }

  template<>
  std::vector<SX> SX::zz_symvar() const {
    SXFunction f("tmp", std::vector<SX>(), make_vector(*this));
    std::vector<SXElement> ret1 = f.free_sx().data();
    std::vector<SX> ret(ret1.size());
    std::copy(ret1.begin(), ret1.end(), ret.begin());
    return ret;
  }

  template<>
  void SX::zz_extractShared(std::vector<SX >& ex,
                            std::vector<SX >& v_sx,
                            std::vector<SX >& vdef_sx,
                            const std::string& v_prefix,
                            const std::string& v_suffix) {

    // Sort the expression
    SXFunction f("tmp", vector<SX>(), ex);

    // Get references to the internal data structures
    const vector<ScalarAtomic>& algorithm = f.algorithm();
    vector<SXElement> work(f.getWorkSize());
    vector<SXElement> work2 = work;

    // Iterator to the binary operations
    vector<SXElement>::const_iterator b_it=f->operations_.begin();

    // Iterator to stack of constants
    vector<SXElement>::const_iterator c_it = f->constants_.begin();

    // Iterator to free variables
    vector<SXElement>::const_iterator p_it = f->free_vars_.begin();

    // Count how many times an expression has been used
    vector<int> usecount(work.size(), 0);

    // Evaluate the algorithm
    vector<SXElement> v, vdef;
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
    for (int i=0; i<vdef.size(); ++i) {
      v_name.str(string());
      v_name << v_prefix << i << v_suffix;
      v.push_back(SXElement::sym(v_name.str()));
    }

    // Mark the above expressions
    for (int i=0; i<vdef.size(); ++i) {
      vdef[i].setTemp(i+1);
    }

    // Save the marked nodes for later cleanup
    vector<SXElement> marked = vdef;

    // Reset iterator
    b_it=f->operations_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_OUTPUT:     ex.at(it->i0).at(it->i2) = work[it->i1];      break;
      case OP_CONST:      work2[it->i0] = work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work2[it->i0] = work[it->i0] = *p_it++; break;
      default:
        {
          switch (it->op) {
            CASADI_MATH_FUN_BUILTIN(work[it->i1], work[it->i2], work[it->i0])
              }
          work2[it->i0] = *b_it++;

          // Replace with intermediate variables
          int ind = work2[it->i0].getTemp()-1;
          if (ind>=0) {
            vdef.at(ind) = work[it->i0];
            work[it->i0] = v.at(ind);
          }
        }
      }
    }

    // Unmark the expressions
    for (vector<SXElement>::iterator it=marked.begin(); it!=marked.end(); ++it) {
      it->setTemp(0);
    }

    // Save v, vdef
    v_sx.resize(v.size());
    std::copy(v.begin(), v.end(), v_sx.begin());
    vdef_sx.resize(vdef.size());
    std::copy(vdef.begin(), vdef.end(), vdef_sx.begin());
  }

  template<>
  SX SX::zz_poly_coeff(const SX& x) const {
    casadi_assert(isscalar());
    casadi_assert(x.isscalar());
    casadi_assert(x.isSymbolic());

    std::vector<SXElement> r;

    SXFunction f("tmp", make_vector(x), make_vector(*this));
    int mult = 1;
    bool success = false;
    for (int i=0; i<1000; ++i) {
      r.push_back((f(SX(0)).at(0)/mult).toScalar());;
      SX j = SX::jac(f);
      if (j.nnz()==0) {
        success = true;
        break;
      }
      f = SXFunction("tmp", make_vector(x), make_vector(j));
      mult*=i+1;
    }

    if (!success) casadi_error("poly: supplied expression does not appear to be polynomial.");

    std::reverse(r.begin(), r.end());

    return r;
  }

  template<>
  SX SX::zz_poly_roots() const {
    const SX& p = *this;
    casadi_assert_message(p.size2()==1,
                          "poly_root(): supplied parameter must be column vector but got "
                          << p.dimString() << ".");
    casadi_assert(p.isdense());
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
      SX ret = vertcat((bm-ds)/a2,
                       (bm+ds)/a2);
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

      SX ret = vertcat(cos(phi/3),
                       cos((phi+2*M_PI)/3),
                       cos((phi+4*M_PI)/3));
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
      SX poly = vertcat(1,
                        f/2,
                        ((f*f -4*h)/16),
                        -g*g/64);
      SX y = poly_roots(poly);

      SX r0 = y(0);
      SX r1 = y(2);

      SX p = sqrt(r0); // two non-zero-roots
      SX q = sqrt(r1);

      SX r = -g/(8*p*q);

      SX s = b/4;

      SX ret = vertcat(p + q + r -s,
                       p - q - r -s,
                       -p + q - r -s,
                       -p - q + r -s);
      return ret;
    } else if (isEqual(p(p.nnz()-1).at(0), 0)) {
      SX ret = vertcat(poly_roots(p(range(p.nnz()-1))),
                       0);
      return ret;
    } else {
      casadi_error("poly_root(): can only solve cases for first or second order polynomial. "
                   "Got order " << p.size1()-1 << ".");
    }

  }

  template<>
  SX SX::zz_eig_symbolic() const {
    const SX& m = *this;
    casadi_assert_message(m.size1()==m.size2(), "eig(): supplied matrix must be square");

    vector<SX> ret;

    /// Bring m in block diagonal form, calculating eigenvalues of each block separately
    std::vector<int> offset;
    std::vector<int> index;
    int nb = m.sparsity().stronglyConnectedComponents(offset, index);

    SX m_perm = m(offset, offset);

    SX l = SX::sym("l");

    for (int k=0; k<nb; ++k) {
      std::vector<int> r = range(index.at(k), index.at(k+1));
      // det(lambda*I-m) = 0
      ret.push_back(poly_roots(poly_coeff(det(SX::eye(r.size())*l-m_perm(r, r)), l)));
    }

    return vertcat(ret);
  }

  SXElement SXElement::zz_simplify() const {
    // Start by expanding the node to a weighted sum
    SX terms, weights;
    expand(*this, weights, terms);

    // Make a scalar product to get the simplified expression
    SX s = mul(terms.T(), weights);
    return s.toScalar();
  }

  template<>
  void SX::printSplit(std::vector<std::string>& nz,
                      std::vector<std::string>& inter) const {
    // Find out which noded can be inlined
    std::map<const SXNode*, int> nodeind;
    for (vector<SXElement>::const_iterator i=begin(); i!=end(); ++i)
      (*i)->can_inline(nodeind);

    // Print expression
    nz.resize(0);
    nz.reserve(nnz());
    inter.resize(0);
    for (vector<SXElement>::const_iterator i=begin(); i!=end(); ++i)
      nz.push_back((*i)->printCompact(nodeind, inter));
  }

  template<> std::vector<SX> SX::get_input(const Function& f) {
    return f.sx_in();
  }

  template<> SX SX::jac(const Function& f, int iind, int oind, bool compact, bool symmetric) {
    return Function(f)->jac_sx(iind, oind, compact, symmetric);
  }

  template<> SX SX::jac(const Function& f, const std::string & iname, int oind,
         bool compact, bool symmetric) {
    return jac(f, f.index_in(iname), oind, compact, symmetric);
  }

  template<> SX SX::jac(const Function& f, int iind, const std::string& oname,
         bool compact, bool symmetric) {
    return jac(f, iind, f.index_out(oname), compact, symmetric);
  }

  template<> SX SX::jac(const Function& f, const std::string& iname, const std::string& oname,
         bool compact, bool symmetric) {
    return jac(f, f.index_in(iname), f.index_out(oname), compact, symmetric);
  }

  template<> SX SX::grad(const Function& f, int iind, int oind) {
    return Function(f)->grad_sx(iind, oind);
  }

  template<> SX SX::grad(const Function& f, const std::string& iname, int oind) {
    return grad(f, f.index_in(iname), oind);
  }

  template<> SX SX::grad(const Function& f, int iind, const std::string& oname) {
    return grad(f, iind, f.index_out(oname));
  }

  template<> SX SX::tang(const Function& f, int iind, int oind) {
    return Function(f)->tang_sx(iind, oind);
  }

  template<> SX SX::grad(const Function& f, const std::string& iname, const std::string& oname) {
    return grad(f, f.index_in(iname), f.index_out(oname));
  }

  template<> SX SX::tang(const Function& f, const std::string& iname, int oind) {
    return tang(f, f.index_in(iname), oind);
  }

  template<> SX SX::tang(const Function& f, int iind, const std::string& oname) {
    return tang(f, iind, f.index_out(oname));
  }

  template<> SX SX::tang(const Function& f, const std::string& iname, const std::string& oname) {
    return tang(f, f.index_in(iname), f.index_out(oname));
  }

  template<> SX SX::hess(const Function& f, int iind, int oind) {
    return Function(f)->hess_sx(iind, oind);
  }

  template<> SX SX::hess(const Function& f, const std::string& iname, int oind) {
    return hess(f, f.index_in(iname), oind);
  }

  template<> SX SX::hess(const Function& f, int iind, const std::string& oname) {
    return hess(f, iind, f.index_out(oname));
  }

  template<> SX SX::hess(const Function& f, const std::string& iname, const std::string& oname) {
    return hess(f, f.index_in(iname), f.index_out(oname));
  }

  template<>
  Function SX::fun(const std::string& name, const Function &f, const Dict& opts) {
    return SXFunction(name, f, opts);
  }

  template<>
  Function SX::fun(const std::string& name, const std::vector<SX>& arg,
                   const std::vector<SX>& res, const Dict& opts) {
    return SXFunction(name, arg, res, opts);
  }

  template<>
  Function SX::fun(const std::string& name,
                   const std::pair< SXDict, std::vector<std::string> >& arg,
                   const std::vector<SX>& res, const Dict& opts) {
    return SXFunction(name, arg, res, opts);
  }

  template<>
  Function SX::fun(const std::string& name, const std::vector<SX>& arg,
                   const std::pair< SXDict, std::vector<std::string> >& res,
                   const Dict& opts) {
    return SXFunction(name, arg, res, opts);
  }

  template<>
  Function SX::fun(const std::string& name,
                   const std::pair< SXDict, std::vector<std::string> >& arg,
                   const std::pair< SXDict, std::vector<std::string> >& res,
                   const Dict& opts) {
    return SXFunction(name, arg, res, opts);
  }

} // namespace casadi
