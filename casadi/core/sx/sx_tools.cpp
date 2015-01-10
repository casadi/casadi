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

  void extractShared(std::vector<SX>& ex, std::vector<SX>& v_sx,
                     std::vector<SX>& vdef_sx,
                     const std::string& v_prefix,
                     const std::string& v_suffix) {

    // Sort the expression
    SXFunction f(vector<SX>(), ex);
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

  void printCompact(const SX& ex, std::ostream &stream) {
    // Extract shared subexpressions from ex
    vector<SX> v, vdef;
    vector<SX> ex_extracted(1, ex);
    extractShared(ex_extracted, v, vdef, "@", "");

    // Print the expression without shared subexpressions
    ex_extracted.at(0).print(stream);

    // Print the shared subexpressions
    if (!v.empty()) {
      stream << endl << "where:" << endl;
      for (int i=0; i<v.size(); ++i) {
        stream << v[i].toScalar() << " := " << vdef[i].toScalar() << endl;
      }
    }
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
    } else if (isEqual(p(p.size()-1).at(0), 0)) {
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
