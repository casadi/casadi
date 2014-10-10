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


#include "sx_tools.hpp"
#include "../function/sx_function_internal.hpp"
#include "../casadi_math.hpp"
#include "../casadi_exception.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../std_vector_tools.hpp"
using namespace std;

namespace casadi {

  SX gauss_quadrature(SX f, const SX &x, const SX &a, const SX &b, int order, const SX& w) {
    casadi_assert_message(order == 5, "gauss_quadrature: order must be 5");
    casadi_assert_message(w.isEmpty(), "gauss_quadrature: empty weights");

    // Change variables to [-1, 1]
    if (!a.toScalar().isEqual(-1) || !b.toScalar().isEqual(1)) {
      SX q1 = (b-a)/2;
      SX q2 = (b+a)/2;

      SXFunction fcn(x, f);
      fcn.init();

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
    SXFunction fcn(x, f);
    vector<SXElement> f_val(5);
    for (int i=0; i<5; ++i)
      f_val[i] = fcn(SX(xi[i])).at(0).toScalar();

    // Weighted sum
    SXElement sum;
    for (int i=0; i<5; ++i)
      sum += wi[i]*f_val[i];

    return sum;
  }

  SX pw_const(const SX &t, const SX &tval, const SX &val) {
    // number of intervals
    int n = val.numel();

    casadi_assert_message(t.isScalar(), "t must be a scalar");
    casadi_assert_message(tval.numel() == n-1, "dimensions do not match");

    SX ret = val.at(0);
    for (int i=0; i<n-1; ++i) {
      ret += (val(0, i+1)-val(0, i)) * (t>=tval(0, i));
    }

    return ret;
  }

  SX pw_lin(const SXElement &t, const SX &tval, const SX &val) {
    // Number of points
    int N = tval.numel();
    casadi_assert_message(N>=2, "pw_lin: N>=2");

    // Gradient for each line segment
    SX g = SX::sparse(1, N-1);
    for (int i=0; i<N-1; ++i) {
      g(0, i) = (val(0, i+1)- val(0, i))/(tval(0, i+1)-tval(0, i));
    }

    // Line segments
    SX lseg = SX::sparse(1, N-1);
    for (int i=0; i<N-1; ++i)
      lseg(0, i) = val(0, i) + g(0, i)*(t-tval(0, i));

    // interior time points
    SX tint = tval(0, range(N-2));

    // Return piecewise linear function
    return pw_const(t, tint, lseg);
  }

  SX if_else(const SX &cond, const SX &if_true, const SX &if_false) {
    return if_else_zero(cond, if_true) + if_else_zero(!cond, if_false);
  }

  SX heaviside(const SX& a) {
    return (1+sign(a))/2;
  }

  SX ramp(const SX& a) {
    return a*heaviside(a);
  }

  SX rectangle(const SX& a) {
    return 0.5*(sign(a+0.5)-sign(a-0.5));
  }

  SX triangle(const SX& a) {
    return rectangle(a.toScalar()/2)*(1-abs(a.toScalar()));
  }

  void simplify(SX &ex) {
    // simplify all non-zero elements
    for (int el=0; el<ex.size(); ++el)
      simplify(ex.at(el));
  }

  void compress(SX &ex, int level) {

    throw CasadiException("SX::compress: Not implemented");

    if (level>0)
      compress(ex, level-1);
  }

  std::vector<SX> substitute(const std::vector<SX> &ex, const std::vector<SX> &v,
                             const std::vector<SX> &vdef) {
    // Assert consistent dimensions
    casadi_assert_warning(v.size()==vdef.size(), "subtitute: number of symbols to replace ( "
                          << v.size() << ") must match number of expressions (" << vdef.size()
                          << ") to replace them with.");

    // Quick return if all equal
    bool all_equal = true;
    for (int k=0; k<v.size(); ++k) {
      if (!v[k].isEqual(vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Check sparsities
    for (int k=0; k<v.size(); ++k) {
      if (v[k].sparsity()!=vdef[k].sparsity()) {
        // Expand vdef to sparsity of v if vdef is scalar
        if (vdef[k].isScalar() && vdef[k].size()==1) {
          std::vector<SX> vdef_mod = vdef;
          vdef_mod[k] = SX(v[k].sparsity(), vdef[k].at(0));
          return substitute(ex, v, vdef_mod);
        } else {
          casadi_error("subsitute(ex, v, vdef): sparsities of v and vdef must match. Got v: "
                       << v[k].dimString() << " and " << "vdef: " << vdef[k].dimString() << ".");
        }
      }
    }


    // Otherwise, evaluate symbolically
    SXFunction F(v, ex);
    F.init();
    return F(vdef);
  }

  SX substitute(const SX &ex, const SX &v, const SX &vdef) {
    return substitute(vector<SX>(1, ex), vector<SX>(1, v), vector<SX>(1, vdef)).front();
  }

  Matrix<double> evalf(const SX &ex, const SX &v, const Matrix<double> &vdef) {
    SXFunction fcn(v, ex);
    fcn.init();
    fcn.input(0).set(vdef);
    fcn.evaluate();
    return fcn.output();
  }

  Matrix<double> evalf(const SX &ex) {
    SXFunction fcn(std::vector< SX >(0), ex);
    fcn.init();
    fcn.evaluate();
    return fcn.output();
  }

  void substituteInPlace(const SX &v, SX &vdef, bool reverse) {
    // Empty vector
    vector<SX> ex;
    substituteInPlace(v, vdef, ex, reverse);
  }

  void substituteInPlace(const SX &v, SX &vdef, std::vector<SX>& ex, bool reverse) {
    casadi_assert_message(v.isSymbolic(), "the variable is not symbolic");
    casadi_assert_message(v.sparsity() == vdef.sparsity(), "the sparsity patterns of the "
                          "expression and its defining bexpression do not match");
    if (v.isEmpty()) return; // quick return if nothing to replace

    // Function inputs
    std::vector<SX> f_in;
    if (!reverse) f_in.push_back(v);

    // Function outputs
    std::vector<SX> f_out;
    f_out.push_back(vdef);
    f_out.insert(f_out.end(), ex.begin(), ex.end());

    // Write the mapping function
    SXFunction f(f_in, f_out);
    f.init();

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
        work[it->i0] = vdef.at(it->i2);
        break;
      case OP_OUTPUT:
        if (it->i0==0) {
          vdef.at(it->i2) = work[it->i1];
          if (reverse) {
            // Use the new variable henceforth, substitute in
            work[it->i1] = v.at(it->i2);
          }
        } else {
          // Auxillary output
          ex[it->i0-1].at(it->i2) = work[it->i1];
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

#if 0
  void replaceDerivatives(SX &ex, const SX &var, const SX &dvar) {
    // Initialize with an empty expression
    SXFunction fcn(ex);

    // Map from the var-node to the new der-node
    std::map<int, SX> dermap;
    for (int i=0; i<var.size(); ++i)
      dermap[fcn.treemap[var[i].get()]] = dvar[i];

    // Replacement map
    std::map<int, SX> replace;

    // Go through all nodes and check if any node is a derivative
    for (int i=0; i<fcn.algorithm.size(); ++i) {
      if (fcn.algorithm[i].op == DER) {

        // find the corresponding derivative
        std::map<int, SX>::iterator r = dermap.find(fcn.algorithm[i].ch0);
        casadi_assert(r != dermap.end());

        replace[i] = r->second;
      }
    }
    SX res;
    SX repres;
    fcn.eval_symbolic(SX(), res, replace, repres);
    ex = res;

    casadi_assert(0);

  }
#endif

#if 0
  void makeSmooth(SX &ex, SX &bvar, SX &bexpr) {
    // Initialize
    SXFunction fcn(SX(), ex);

    casadi_assert(bexpr.isEmpty());

    // Nodes to be replaced
    std::map<int, SXElement> replace;

    // Go through all nodes and check if any node is non-smooth
    for (int i=0; i<fcn->algorithm.size(); ++i) {

      // Check if we have a step node
      if (fcn->algorithm[i].op == STEP) {

        // Get the index of the child
        int ch0 = fcn->algorithm[i].ch[0];

        // Binary variable corresponding to the the switch
        SX sw;

#if 0
        // Find out if the switch has already been added
        for (int j=0; j<bexpr.size(); ++j)
          if (bexpr[j].isEqual(algorithm[i]->child0)) {
            sw = bvar[j];
            break;
          }
#endif

        if (sw.isEmpty()) { // the switch has not yet been added
          // Get an approriate name of the switch
          std::stringstream name;
          name << "sw_" << bvar.size2();
          sw = SXElement::sym(name.str());

          // Add to list of switches
          bvar << sw;
          //        bexpr << algorithm[i]->child0;
        }

        // Add to the substition map
        replace[i] = sw[0];
      }
    }
    SX res;
    fcn->eval(SX(), res, replace, bexpr);

    for (int i=0; i<bexpr.size(); ++i)
      bexpr[i] = bexpr[i]->dep(0);

    ex = res;

#if 0
    // Make sure that the binding expression is smooth
    bexpr.init(SX());
    SX b;
    bexpr.eval_symbolic(SX(), b, replace, bexpr);
    bexpr = b;
#endif
  }
#endif

  SX spy(const SX& A) {
    SX s = SX::sparse(A.size1(), A.size2());
    for (int i=0; i<A.size2(); ++i)
      for (int j=0; j<A.size1(); ++j)
        if (!A(j, i).toScalar()->isZero())
          s(j, i) = 1;
    return s;
  }

  bool dependsOn(const SX& ex, const SX &arg) {
    if (ex.size()==0) return false;

    // Construct a temporary algorithm
    SXFunction temp(arg, ex);
    temp.init();
    temp.spInit(true);

    bvec_t* input_ =  get_bvec_t(temp.input().data());
    // Make a column with all variables active
    std::fill(input_, input_+temp.input().size(), bvec_t(1));
    bvec_t* output_ = get_bvec_t(temp.output().data());
    // Perform a single dependency sweep
    temp.spEvaluate(true);

    // Loop over results
    for (int i=0;i<temp.output().size();++i) {
      if (output_[i]) return true;
    }

    return false;
  }

  SX gradient(const SX& ex, const SX &arg) {
    SXFunction temp(arg, ex); // make a runtime
    temp.init();
    return temp.grad();
  }

  SX tangent(const SX& ex, const SX &arg) {
    SXFunction temp(arg, ex); // make a runtime
    temp.init();
    return temp.tang();
  }

  SX jacobian(const SX& ex, const SX &arg) {
    SXFunction temp(arg, ex); // make a runtime
    temp.init();
    return temp.jac();
  }

  void hessian(const SX& ex, const SX &arg, SX &H, SX &g) {
    g = gradient(ex, arg);

    SXFunction temp(arg, g); // make a runtime
    temp.init();
    H = temp.jac(0, 0, false, true);
  }

  SX hessian(const SX& ex, const SX &arg) {
    SX H, g;
    hessian(ex, arg, H, g);
    return H;
  }

  void expand(const SX& ex2, SX &ww, SX& tt) {
    casadi_assert(ex2.isScalar());
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

  void simplify(SXElement& ex) {
    // Start by expanding the node to a weighted sum
    SX terms, weights;
    expand(ex, weights, terms);

    // Make a scalar product to get the simplified expression
    SX s = mul(terms.T(), weights);
    ex = s.toScalar();
  }

  SX taylor(const SX& ex, const SX& x, const SX& a, int order) {
    casadi_assert(x.isScalar() && a.isScalar());
    if (ex.size()!=ex.numel())
      throw CasadiException("taylor: not implemented for sparse matrices");
    SX ff = vec(ex.T());

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
    return reshape(result, ex.size2(), ex.size1()).T();
  }

  SX mtaylor(const SX& ex, const SX& x, const SX& around, int order) {
    return mtaylor(ex, x, around, order, std::vector<int>(x.size(), 1));
  }

  /// \cond
  SX mtaylor_recursive(const SX& ex, const SX& x, const SX& a, int order,
                       const std::vector<int>&order_contributions,
                       const SXElement & current_dx=casadi_limits<SXElement>::one,
                       double current_denom=1, int current_order=1) {
    SX result = substitute(ex, x, a)*current_dx/current_denom;
    for (int i=0;i<x.size();i++) {
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
  /// \endcond

  SX mtaylor(const SX& ex, const SX& x, const SX& a, int order,
             const std::vector<int>&order_contributions) {
    casadi_assert_message(ex.size()==ex.numel() && x.size()==x.numel(),
                          "mtaylor: not implemented for sparse matrices");

    casadi_assert_message(x.size()==order_contributions.size(),
                          "mtaylor: number of non-zero elements in x (" <<  x.size()
                          << ") must match size of order_contributions ("
                          << order_contributions.size() << ")");

    return reshape(
             mtaylor_recursive(vec(ex), x, a, order, order_contributions),
             ex.size2(), ex.size1()).T();
  }

  int countNodes(const SX& A) {
    SXFunction f(SX(), A);
    f.init();
    return f.countNodes();
  }


  std::string getOperatorRepresentation(const SXElement& x, const std::vector<std::string>& args) {
    if (!x.hasDep())
        throw CasadiException("getOperatorRepresentation: SXElement must be binary operator");
    if (args.size() == 0 || (casadi_math<double>::ndeps(x.getOp())==2 && args.size() < 2))
        throw CasadiException("getOperatorRepresentation: not enough arguments supplied");
    std::stringstream s;
    casadi_math<double>::print(x.getOp(), s, args[0], args[1]);
    return s.str();
  }

  void makeSemiExplicit(const SX& f, const SX& x, SX& fe, SX& fi, SX& xe, SX& xi) {
    casadi_assert(f.isDense());
    casadi_assert(x.isDense());

    // Create the implicit function
    SXFunction fcn(x, f);
    fcn.init();

    // Get the sparsity pattern of the Jacobian (no need to actually form the Jacobian)
    Sparsity Jsp = fcn.jacSparsity();

    // Free the function
    fcn = SXFunction();

    // Make a BLT sorting of the Jacobian (a Dulmage-Mendelsohn decomposition)
    std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    Jsp.dulmageMendelsohn(colperm, rowperm, colperm, rowblock, coarse_colblock, coarse_rowblock);

    // Make sure that the Jacobian is full rank
    casadi_assert(coarse_colblock[0]==0);
    casadi_assert(coarse_colblock[1]==0);
    casadi_assert(coarse_colblock[2]==0);
    casadi_assert(coarse_colblock[3]==coarse_colblock[4]);

    casadi_assert(coarse_rowblock[0]==0);
    casadi_assert(coarse_rowblock[1]==0);
    casadi_assert(coarse_rowblock[2]==coarse_rowblock[3]);
    casadi_assert(coarse_rowblock[3]==coarse_rowblock[4]);

    // Permuted equations
    vector<SXElement> fp(f.size());
    for (int i=0; i<fp.size(); ++i) {
      fp[i] = f.elem(0, colperm[i]);
    }

    // Permuted variables
    vector<SXElement> xp(x.size());
    for (int i=0; i<xp.size(); ++i) {
      xp[i]= x.elem(0, rowperm[i]);
    }

    // Number of blocks
    int nb = colblock.size()-1;

    // Block equations
    vector<SXElement> fb;

    // Block variables
    vector<SXElement> xb;

    // Block variables that enter linearly and nonlinearily respectively
    vector<SXElement> xb_lin, xb_nonlin;

    // The separated variables and equations
    vector<SXElement> fev, fiv, xev, xiv;

    // Loop over blocks
    for (int b=0; b<nb; ++b) {

      // Get the local equations
      fb.clear();
      for (int i=colblock[b]; i<colblock[b+1]; ++i) {
        fb.push_back(fp[i]);
      }

      // Get the local variables
      xb.clear();
      for (int i=rowblock[b]; i<rowblock[b+1]; ++i) {
        xb.push_back(xp[i]);
      }

      // We shall find out which variables enter nonlinearily in the equations,
      // for this we need a function that will depend on all the variables
      SXFunction fcnb_all(xb, inner_prod(SX(fb), SX::sym("dum1", fb.size())));
      fcnb_all.init();

      // Take the gradient of this function to find out
      // which variables enter in the function (should be all)
      SX fcnb_dep = fcnb_all.grad();

      // Make sure that this expression is dense (otherwise, some variables would not enter)
      casadi_assert(fcnb_dep.isDense());

      // Multiply this expression with a new dummy vector and take the jacobian to
      // find out which variables enter nonlinearily
      SXFunction fcnb_nonlin(xb, inner_prod(fcnb_dep, SX::sym("dum2", fcnb_dep.size())));
      fcnb_nonlin.init();
      Sparsity sp_nonlin = fcnb_nonlin.jacSparsity().transpose();

      // Get the subsets of variables that appear nonlinearily
      vector<bool> nonlin(sp_nonlin.size1(), false);
      for (int el=0; el<sp_nonlin.size(); ++el) {
        nonlin[sp_nonlin.row(el)] = true;
      }
      /*    cout << "nonlin = " << nonlin << endl;*/

      // Separate variables
      xb_lin.clear();
      xb_nonlin.clear();
      for (int i=0; i<nonlin.size(); ++i) {
        if (nonlin[i])
          xb_nonlin.push_back(xb[i]);
        else
          xb_lin.push_back(xb[i]);
      }

      // If there are only nonlinear variables
      if (xb_lin.empty()) {
        // Substitute the already determined variables
        fb = substitute(SX(fb), SX(xev), SX(fev)).data();

        // Add to the implicit variables and equations
        fiv.insert(fiv.end(), fb.begin(), fb.end());
        xiv.insert(xiv.end(), xb.begin(), xb.end());
      } else {
        // Write the equations as a function of the linear variables
        SXFunction fcnb(xb_lin, fb);
        fcnb.init();

        // Write the equation in matrix form
        SX Jb = fcnb.jac();
        SX rb = -fcnb(SX::zeros(1, xb_lin.size())).at(0);

        // Simple solve if there are no nonlinear variables
        if (xb_nonlin.empty()) {

          // Check if 1-by-1 block
          if (Jb.numel()==1) {
            // Simple division if Jb scalar
            rb /= Jb;
          } else {
            // Solve system of equations
            rb = solve(Jb, rb.T()).T();
          }

          // Substitute the already determined variables
          rb = substitute(rb, SX(xev), SX(fev));

          // Add to the explicit variables and equations
          fev.insert(fev.end(), rb.begin(), rb.end());
          xev.insert(xev.end(), xb.begin(), xb.end());

        } else { // There are both linear and nonlinear variables

          // Make a Dulmage-Mendelsohn decomposition
          std::vector<int> rowpermb, colpermb, rowblockb, colblockb,
              coarse_rowblockb, coarse_colblockb;
          Jb.sparsity().dulmageMendelsohn(rowpermb, colpermb, rowblockb, colblockb,
                                          coarse_rowblockb, coarse_colblockb);

          Matrix<int>(Jb.sparsity(), 1).printDense();
          Jb.printDense();






          cout << colpermb << endl;
          cout << rowpermb << endl;
          cout << colblockb << endl;
          cout << rowblockb << endl;
          cout << coarse_colblockb << endl;
          cout << coarse_rowblockb << endl;

          casadi_warning("tearing not implemented");


          // Substitute the already determined variables
          fb = substitute(SX(fb), SX(xev), SX(fev)).data();

          // Add to the implicit variables and equations
          fiv.insert(fiv.end(), fb.begin(), fb.end());
          xiv.insert(xiv.end(), xb.begin(), xb.end());

        }
      }
    }

    fi = SX(fiv);
    fe = SX(fev);
    xi = SX(xiv);
    xe = SX(xev);
  }

  SX getFree(const SX& ex) {
    SXFunction f(vector<SX>(), ex);
    f.init();
    return f.getFree();
  }

  SX jacobianTimesVector(const SX &ex, const SX &arg, const SX &v, bool transpose_jacobian) {
    SXFunction f(arg, ex);
    f.init();

    // Split up v
    vector<SX> vv = horzsplit(v);

    // Make sure well-posed
    casadi_assert(vv.size() >= 1);
    casadi_assert(ex.isVector());
    casadi_assert(arg.isVector());
    if (transpose_jacobian) {
      casadi_assert(v.size1()==ex.size1());
    } else {
      casadi_assert(v.size1()==arg.size1());
    }

    // Number of sensitivities
    int nfsens = transpose_jacobian ? 0 : vv.size();
    int nasens = transpose_jacobian ? vv.size() : 0;

    // Assemble arguments and directional derivatives
    vector<SX> argv = f.inputExpr();
    vector<SX> resv = f.outputExpr();
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

  void extractShared(std::vector<SXElement>& ex, std::vector<SXElement>& v,
                     std::vector<SXElement>& vdef,
                     const std::string& v_prefix, const std::string& v_suffix) {

    // Sort the expression
    SXFunction f(vector<SX>(), vector<SX>(1, ex));
    f.init();

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

    // New variables and definitions
    v.clear();
    vdef.clear();

    // Evaluate the algorithm
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
      case OP_OUTPUT:     ex[it->i2] = work[it->i1];      break;
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
  }

  void printCompact(const SX& ex, std::ostream &stream) {
    // Extract shared subexpressions from ex
    vector<SXElement> v, vdef;
    SX ex_extracted = ex;
    extractShared(ex_extracted.data(), v, vdef, "@", "");

    // Print the expression without shared subexpressions
    ex_extracted.print(stream);

    // Print the shared subexpressions
    if (!v.empty()) {
      stream << endl << "where:" << endl;
      for (int i=0; i<v.size(); ++i) {
        stream << v[i] << " := " << vdef[i] << endl;
      }
    }
  }

  void substituteInPlace(const std::vector<SX>& v, std::vector<SX>& vdef, std::vector<SX>& ex,
                         bool reverse) {
    casadi_assert(v.size()==vdef.size());

    // Quick return if empty or single expression
    if (v.empty()) {
      return;
    } else if (v.size()==1) {
      substituteInPlace(v.front(), vdef.front(), ex, reverse);
      return;
    }

    // Count number of scalar variables
    int n =0;
    for (int i=0; i<v.size(); ++i) {
      casadi_assert_message(v[i].sparsity() == vdef[i].sparsity(),
                            "the sparsity patterns of the expression and its "
                            "defining expression do not match");
      n += v[i].size();
    }

    // Gather all variables
    SX v_all = SX::zeros(1, n);
    SX vdef_all = SX::zeros(1, n);
    vector<SXElement>::iterator it_v = v_all.begin();
    vector<SXElement>::iterator it_vdef = vdef_all.begin();
    for (int i=0; i<v.size(); ++i) {
      int nv = v[i].size();
      copy(v[i].begin(), v[i].end(), it_v);
      copy(vdef[i].begin(), vdef[i].end(), it_vdef);
      it_v += nv;  it_vdef += nv;
    }

    // Substitute
    substituteInPlace(v_all, vdef_all, ex, reverse);

    // Collect the result
    it_vdef = vdef_all.begin();
    for (int i=0; i<v.size(); ++i) {
      int nv = v[i].size();
      copy(it_vdef, it_vdef+nv, vdef[i].begin());
      it_vdef += nv;
    }
  }

  std::vector<SXElement> getSymbols(const SX& e) {
    SXFunction f(std::vector<SX>(), e);
    f.init();
    return f.getFree();
  }

  SX poly_coeff(const SX& ex, const SX&x) {
    casadi_assert(ex.isScalar());
    casadi_assert(x.isScalar());
    casadi_assert(x.isSymbolic());

    SX ret;

    SXFunction f(x, ex);
    f.init();
    int mult = 1;
    bool success = false;
    for (int i=0;i<1000;++i) {
      ret.append(f(SX(0)).at(0)/mult);
      SX j = f.jac();
      if (j.size()==0) {
        success = true;
        break;
      }
      f = SXFunction(x, j);
      f.init();
      mult*=i+1;
    }

    if (!success) casadi_error("poly: supplied expression does not appear to be polynomial.");

    std::reverse(ret.data().begin(), ret.data().end());

    return ret;


  }

  SX poly_roots(const SX& p) {
    casadi_assert_message(p.size2()==1,
                          "poly_root(): supplied parameter must be column vector but got "
                          << p.dimString() << ".");
    casadi_assert(p.isDense());
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
      SX ret;
      ret.append((bm-ds)/a2);
      ret.append((bm+ds)/a2);
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

      SX ret;
      ret.append(cos(phi/3));
      ret.append(cos((phi+2*M_PI)/3));
      ret.append(cos((phi+4*M_PI)/3));
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
      SX poly;
      poly.append(1);
      poly.append(f/2);
      poly.append((f*f -4*h)/16);
      poly.append(-g*g/64);
      SX y = poly_roots(poly);

      SX r0 = y(0);
      SX r1 = y(2);

      SX p = sqrt(r0); // two non-zero-roots
      SX q = sqrt(r1);

      SX r = -g/(8*p*q);

      SX s = b/4;

      SX ret;
      ret.append(p + q + r -s);
      ret.append(p - q - r -s);
      ret.append(-p + q - r -s);
      ret.append(-p - q + r -s);

      return ret;
    } else if (p(p.size()-1).at(0).isEqual(0)) {
      SX ret = poly_roots(p(range(p.size()-1)));
      ret.append(0);
      return ret;
    } else {
      casadi_error("poly_root(): can only solve cases for first or second order polynomial. "
                   "Got order " << p.size1()-1 << ".");
    }

  }

  SX eig_symbolic(const SX& m) {
    casadi_assert_message(m.size1()==m.size2(), "eig(): supplied matrix must be square");

    SX ret;

    /// Bring m in block diagonal form, calculating eigenvalues of each block separately
    std::vector<int> offset;
    std::vector<int> index;
    int nb = m.sparsity().stronglyConnectedComponents(offset, index);

    SX m_perm = m(offset, offset);

    SX l = SX::sym("l");

    for (int k=0;k<nb;++k) {
      std::vector<int> r = range(index.at(k), index.at(k+1));
      // det(lambda*I-m) = 0
      ret.append(poly_roots(poly_coeff(det(SX::eye(r.size())*l-m_perm(r, r)), l)));
    }

    return ret;
  }

} // namespace casadi
