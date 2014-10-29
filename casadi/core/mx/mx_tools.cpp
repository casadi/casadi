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


#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/mx_function.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../std_vector_tools.hpp"
#include "../function/mx_function_internal.hpp"

using namespace std;

namespace casadi {

  // Helper function
  bool has_empty(const vector<MX>& x, bool both=false) {
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (i->isEmpty(both)) return true;
    }
    return false;
  }

  vector<MX> trim_empty(const vector<MX>& x, bool both=false) {
    vector<MX> ret;
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (!i->isEmpty(both)) ret.push_back(*i);
    }
    return ret;
  }

  MX horzcat(const vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      vector<MX> ne = trim_empty(x, true);
      for (int i=0;i<ne.size();i++) {
        casadi_assert_message(ne[i].size1()==ne[0].size1(),
                      "horzcat dimension mismatch  " <<
                      "x[" << i << "]:" << ne[i].dimString() <<
                      " and x[0]: " << ne[0].dimString() << ".");
      }
    }

    if (x.empty()) {
      return MX();
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want horzcat(zeros(0,5),zeros(0,5)) -> zeros(0,10)
        ret = trim_empty(x, true);
        int s = 0;
        for (int i=0;i<ret.size();++i) {
          s+= ret[i].size2();
        }
        return MX::zeros(0, s);
      } else {
        return horzcat(ret);
      }
    } else {
      return x.front()->getHorzcat(x);
    }
  }

  MX diagcat(const vector<MX>& x) {
    if (x.empty()) {
      return MX();
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want diagcat(zeros(5,0),zeros(5,0)) -> zeros(10,0)
        ret = trim_empty(x, true);
        int s1 = 0;
        int s2 = 0;
        for (int i=0;i<ret.size();++i) {
          s1+= ret[i].size1();
          s2+= ret[i].size2();
        }
        return MX::zeros(s1, s2);
      } else {
        return diagcat(ret);
      }
    } else {
      return x.front()->getDiagcat(x);
    }
  }

  MX vertcat(const vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      vector<MX> ne = trim_empty(x, true);
      for (int i=0;i<ne.size();i++) {
        casadi_assert_message(ne[i].size2()==ne[0].size2(),
                      "vertcat dimension mismatch  " <<
                      "x[" << i << "]:" << ne[i].dimString() <<
                      " and x[0]: " << ne[0].dimString() << ".");
      }
    }

    if (x.empty()) {
      return MX();
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want vertcat(zeros(5,0),zeros(5,0)) -> zeros(10,0)
        ret = trim_empty(x, true);
        int s = 0;
        for (int i=0;i<ret.size();++i) {
          s+= ret[i].size1();
        }
        return MX::zeros(s, 0);
      } else {
        return vertcat(ret);
      }
    } else if (!x.front().isVector()) {
      // Vertcat operation only supports vectors, rewrite using horzcat
      vector<MX> xT = x;
      for (vector<MX>::iterator i=xT.begin(); i!=xT.end(); ++i) *i = i->T();
      return horzcat(xT).T();
    } else {
      return x.front()->getVertcat(x);
    }
  }

  std::vector<MX> horzsplit(const MX& x, const std::vector<int>& offset) {
    // Consistency check
    casadi_assert(offset.size()>=1);
    casadi_assert(offset.front()==0);
    casadi_assert(offset.back()==x.size2());
    casadi_assert(isMonotone(offset));

    // Trivial return if possible
    if (offset.size()==1) {
      return vector<MX>(0);
    } else if (offset.size()==2) {
      return vector<MX>(1, x);
    } else {
      return x->getHorzsplit(offset);
    }
  }

  std::vector<MX> horzsplit(const MX& x, int incr) {
    casadi_assert(incr>=1);
    vector<int> offset2 = range(0, x.size2(), incr);
    offset2.push_back(x.size2());
    return horzsplit(x, offset2);
  }

  std::vector<MX> diagsplitNative(const MX& x,
    const std::vector<int>& offset1, const std::vector<int>& offset2) {
    // Consistency check
    casadi_assert(offset1.size()>=1);
    casadi_assert(offset1.front()==0);
    casadi_assert(offset1.back()==x.size1());
    casadi_assert(isMonotone(offset1));

    // Consistency check
    casadi_assert(offset2.size()>=1);
    casadi_assert(offset2.front()==0);
    casadi_assert(offset2.back()==x.size2());
    casadi_assert(isMonotone(offset2));

    return x->getDiagsplit(offset1, offset2);
  }

  std::vector<MX> vertsplit(const MX& x, const std::vector<int>& offset) {
    if (x.isVector()) {
      // Consistency check
      casadi_assert(offset.size()>=1);
      casadi_assert(offset.front()==0);
      casadi_assert(offset.back()==x.size1());
      casadi_assert(isMonotone(offset));

      // Trivial return if possible
      if (offset.size()==1) {
        return vector<MX>();
      } else if (offset.size()==2) {
        return vector<MX>(1, x);
      } else {
        return x->getVertsplit(offset);
      }
    } else {
      std::vector<MX> ret = horzsplit(x.T(), offset);
      MX (*transposeMX)(const MX& x) = transpose;
      std::transform(ret.begin(), ret.end(), ret.begin(), transposeMX);
      return ret;
    }
  }

  std::vector<MX> vertsplit(const MX& x, int incr) {
    casadi_assert(incr>=1);
    vector<int> offset1 = range(0, x.size1(), incr);
    offset1.push_back(x.size1());
    return vertsplit(x, offset1);
  }

  std::vector< std::vector<MX> > blocksplit(const MX& x, const std::vector<int>& vert_offset,
                                            const std::vector<int>& horz_offset) {
    std::vector<MX> rows = vertsplit(x, vert_offset);
    std::vector< std::vector<MX> > ret;
    for (int i=0;i<rows.size();++i) {
      ret.push_back(horzsplit(rows[i], horz_offset));
    }
    return ret;
  }

  std::vector< std::vector<MX > > blocksplit(const MX& x, int vert_incr, int horz_incr) {
    casadi_assert(horz_incr>=1);
    casadi_assert(vert_incr>=1);
    vector<int> offset1 = range(0, x.size1(), vert_incr);
    offset1.push_back(x.size1());
    vector<int> offset2 = range(0, x.size2(), horz_incr);
    offset2.push_back(x.size2());
    return blocksplit(x, offset1, offset2);
  }

  MX horzcat(const MX& a, const MX& b) {
    vector<MX> ab;
    ab.push_back(a);
    ab.push_back(b);
    return horzcat(ab);
  }

  MX vertcat(const MX& a, const MX& b) {
    vector<MX> ab;
    ab.push_back(a);
    ab.push_back(b);
    return vertcat(ab);
  }

  MX veccat(const vector<MX>& comp) {
    MX (&f)(const MX&) = vec;
    return vertcat(applymap(f, comp));
  }

  MX vecNZcat(const vector<MX>& comp) {
    return vertcat(applymap(vecNZ, comp));
  }

  MX norm_2(const MX &x) {
    if (x.isVector()) {
      return norm_F(x);
    } else {
      return x->getNorm2();
    }
  }

  MX norm_F(const MX &x) {
    return x->getNormF();
  }

  MX norm_1(const MX &x) {
    return x->getNorm1();
  }

  MX norm_inf(const MX &x) {
    return x->getNormInf();
  }

  MX mul(const MX &x, const MX &y, const Sparsity& sp_z) {
    return x.mul(y, sp_z);
  }

  MX mul(const std::vector< MX > &args) {
    casadi_assert_message(args.size()>=1,
                          "mul(std::vector< MX > &args): supplied list must not be empty.");
    if (args.size()==1) return args[0];
    MX ret = args[0].mul(args[1]);
    for (int i=2;i<args.size();++i) {
      ret = ret.mul(args[i]);
    }
    return ret;
  }

  MX inner_prod(const MX &x, const MX &y) {
    return x.inner_prod(y);
  }

  MX outer_prod(const MX &x, const MX &y) {
    return x.outer_prod(y);
  }

  void simplify(MX& ex) {
    if (!ex.isEmpty(true)) {
      ex->simplifyMe(ex);
    }
  }

  MX transpose(const MX &x) {
    return x.T();
  }

  MX reshape(const MX &x, std::pair<int, int> rc) {
    return reshape(x, rc.first, rc.second);
  }

  MX reshape(const MX &x, int nrow, int ncol) {
    if (nrow==x.size1() && ncol==x.size2())
      return x;
    else
      return reshape(x, x.sparsity().reshape(nrow, ncol));
  }

  MX reshape(const MX &x, const Sparsity& sp) {
    // quick return if already the right shape
    if (sp==x.sparsity())
      return x;

    // make sure that the patterns match
    casadi_assert(sp.isReshape(x.sparsity()));

    // Create a reshape node
    return x->getReshape(sp);
  }

  MX vec(const MX& x) {
    if (x.isVector()) {
      return x;
    } else {
      return reshape(x, x.numel(), 1);
    }
  }

  MX vecNZ(const MX& x) {
    if (x.isDense()) {
      return vec(x);
    } else {
      return x->getGetNonzeros(Sparsity::dense(x.size(), 1), range(x.size()));
    }
  }

  MX if_else(const MX &cond, const MX &if_true, const MX &if_false) {
    return if_else_zero(cond, if_true) + if_else_zero(!cond, if_false);
  }

  MX unite(const MX& A, const MX& B) {
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    Sparsity sp = A.sparsity().patternUnion(B.sparsity(), mapping);

    // Split up the mapping
    std::vector<int> nzA, nzB;

    // Copy sparsity
    for (int k=0; k<mapping.size(); ++k) {
      if (mapping[k]==1) {
        nzA.push_back(k);
      } else if (mapping[k]==2) {
        nzB.push_back(k);
      } else {
        throw CasadiException("Pattern intersection not empty");
      }
    }

    // Create mapping
    MX ret = MX::zeros(sp);
    ret = A->getSetNonzeros(ret, nzA);
    ret = B->getSetNonzeros(ret, nzB);
    return ret;
  }

  MX trace(const MX& A) {
    casadi_assert_message(A.size2() == A.size1(), "trace: must be square");
    MX res(0);
    for (int i=0; i < A.size2(); i ++) {
      res+=A(i, i);
    }
    return res;
  }

  MX repmat(const MX &A, int n, int m) {
    // Quick return if possible
    if (n==1 &&  m==1)
      return A;

    // First concatenate horizontally
    MX col = horzcat(std::vector<MX >(m, A));

    // Then vertically
    return vertcat(std::vector<MX >(n, col));
  }

  /**
     MX clip(const MX& A, const Sparsity& sp) {
     // Join the sparsity patterns
     std::vector<int> mapping;
     Sparsity sp = A.sparsity().patternIntersection(sp, mapping);

     // Split up the mapping
     std::vector<int> nzA, nzB;

     // Copy sparsity
     for (int k=0; k<mapping.size(); ++k) {
     if (mapping[k]<0) {
     nzA.push_back(k);
     } else if (mapping[k]>0) {
     nzB.push_back(k);
     } else {
     throw CasadiException("Pattern intersection not empty");
     }
     }

     // Create mapping
     MX ret;
     ret.assignNode(new Mapping(sp));
     ret->assign(A, range(nzA.size()), nzA);
     ret->assign(B, range(nzB.size()), nzB);
     return ret;

     }
  */

  MX dense(const MX& x) {
    MX ret = x;
    ret.densify();
    return ret;
  }

  MX createParent(std::vector<MX> &deps) {
    // First check if arguments are symbolic
    for (int k=0;k<deps.size();k++) {
      if (!deps[k].isSymbolic())
          throw CasadiException("createParent: the argumenst must be pure symbolic");
    }

    // Collect the sizes of the depenencies
    std::vector<int> index(deps.size()+1, 0);
    for (int k=0;k<deps.size();k++) {
      index[k+1] =  index[k] + deps[k].size();
    }

    // Create the parent
    MX P = MX::sym("P", index[deps.size()], 1);

    std::vector<MX> Ps = vertsplit(P, index);

    // Make the arguments dependent on the parent
    for (int k=0;k<deps.size();k++) {
      deps[k] = MX(deps[k].sparsity(), Ps[k]);
    }

    return P;
  }

  MX createParent(const std::vector<Sparsity> &deps, std::vector<MX>& children) {
    // Collect the sizes of the depenencies
    std::vector<int> index(deps.size()+1, 0);
    for (int k=0;k<deps.size();k++) {
      index[k+1] =  index[k] + deps[k].size();
    }

    // Create the parent
    MX P = MX::sym("P", index[deps.size()], 1);
    std::vector<MX> Ps = vertsplit(P, index);

    children.resize(deps.size());

    // Make the arguments dependent on the parent
    for (int k=0;k<deps.size();k++) {
      children[k] =  MX(deps[k], Ps[k]);
    }

    return P;
  }

  MX createParent(const std::vector<MX> &deps, std::vector<MX>& children) {
    children = deps;
    MX P = createParent(children);
    return P;
  }

  MX diag(const MX& x) {
    // Nonzero mapping
    std::vector<int> mapping;

    // Get the sparsity
    Sparsity sp = x.sparsity().getDiag(mapping);

    // Create a reference to the nonzeros
    return x->getGetNonzeros(sp, mapping);
  }

  MX blkdiag(const std::vector<MX> &A) {
    return diagcat(A);
  }

  MX blkdiag(const MX &A, const MX& B) {
    std::vector<MX> ret;
    ret.push_back(A);
    ret.push_back(B);
    return blkdiag(ret);
  }

  int countNodes(const MX& A) {
    MXFunction f(vector<MX>(), A);
    f.init();
    return f.countNodes();
  }

  MX sumCols(const MX &x) {
    return mul(x, MX::ones(x.size2(), 1));
  }

  MX sumRows(const MX &x) {
    return mul(MX::ones(1, x.size1()), x);
  }

  MX sumAll(const MX &x) {
    return sumRows(sumCols(x));
  }


  MX polyval(const MX& p, const MX& x) {
    casadi_assert_message(p.isDense(), "polynomial coefficients vector must be a vector");
    casadi_assert_message(p.isVector() && p.size()>0, "polynomial coefficients must be a vector");
    MX ret = p[0];
    for (int i=1; i<p.size(); ++i) {
      ret = ret*x + p[i];
    }
    return ret;
  }

  std::string getOperatorRepresentation(const MX& x, const std::vector<std::string>& args) {
    std::stringstream s;
    const MXNode* node = dynamic_cast<const MXNode*>(x.get());
    node->printPart(s, 0);
    if (x.isUnary()) {
      s << args[0];
      node->printPart(s, 1);
    } else {
      for (int i=0;i<args.size();++i) {
        s << args[i];
        node->printPart(s, 1+i);
      }
    }

    return s.str();
  }

  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef, bool reverse) {
    // Empty vector
    vector<MX> ex;
    substituteInPlace(v, vdef, ex, reverse);
  }

  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef,
                         std::vector<MX>& ex, bool reverse) {
    casadi_assert_message(v.size()==vdef.size(),
                          "Mismatch in the number of expression to substitute.");
    for (int k=0; k<v.size(); ++k) {
      casadi_assert_message(v[k].isSymbolic(), "Variable " << k << " is not symbolic");
      casadi_assert_message(v[k].sparsity() == vdef[k].sparsity(),
                            "Inconsistent sparsity for variable " << k << ".");
    }
    casadi_assert_message(reverse==false, "Not implemented");

    // quick return if nothing to replace
    if (v.empty()) return;

    // Function inputs
    std::vector<MX> f_in = v;

    // Function outputs
    std::vector<MX> f_out = vdef;
    f_out.insert(f_out.end(), ex.begin(), ex.end());

    // Write the mapping function
    MXFunction f(f_in, f_out);
    f.init();

    // Get references to the internal data structures
    std::vector<MXAlgEl>& algorithm = f->algorithm_;
    vector<MX> work(f.getWorkSize());
    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;

    for (vector<MXAlgEl>::iterator it=algorithm.begin(); it!=algorithm.end(); ++it) {
      switch (it->op) {
      case OP_INPUT:
        work.at(it->res.front()) = vdef.at(it->arg.front());
        break;
      case OP_PARAMETER:
      case OP_CONST:
        work.at(it->res.front()) = it->data;
        break;
      case OP_OUTPUT:
        if (it->res.front()<vdef.size()) {
          vdef.at(it->res.front()) = work.at(it->arg.front());
        } else {
          ex.at(it->res.front()-vdef.size()) = work.at(it->arg.front());
        }
        break;
      default:
        {
          input_p.resize(it->arg.size());
          for (int i=0; i<input_p.size(); ++i) {
            int el = it->arg[i];
            input_p[i] = el<0 ? 0 : &work.at(el);
          }

          output_p.resize(it->res.size());
          for (int i=0; i<output_p.size(); ++i) {
            int el = it->res[i];
            output_p[i] = el<0 ? 0 : &work.at(el);
          }

          it->data->evaluateMX(input_p, output_p, dummy_p, dummy_p, dummy_p, dummy_p, false);
        }
      }
    }
  }

  MX substitute(const MX &ex, const MX& v, const MX& vdef) {
    return substitute(vector<MX>(1, ex), vector<MX>(1, v), vector<MX>(1, vdef)).front();
  }

  std::vector<MX> substitute(const std::vector<MX> &ex, const std::vector<MX> &v,
                             const std::vector<MX> &vdef) {
    // Assert consistent dimensions
    casadi_assert(v.size()==vdef.size());

    // Quick return if all equal
    bool all_equal = true;
    for (int k=0; k<v.size(); ++k) {
      if (!v[k].isEqual(vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Otherwise, evaluate symbolically
    MXFunction F(v, ex);
    F.init();
    return F.call(vdef, true);
  }

  MX graph_substitute(const MX &ex, const std::vector<MX> &v, const std::vector<MX> &vdef) {
    return graph_substitute(std::vector<MX>(1, ex), v, vdef).at(0);
  }

  std::vector<MX> graph_substitute(const std::vector<MX> &ex, const std::vector<MX> &expr,
                                   const std::vector<MX> &exprs) {
    casadi_assert_message(expr.size()==exprs.size(),
                          "Mismatch in the number of expression to substitute: "
                          << expr.size() << " <-> " << exprs.size() << ".");

    // Sort the expression
    MXFunction f(vector<MX>(), ex);
    f.init();

    // Get references to the internal data structures
    const vector<MXAlgEl>& algorithm = f.algorithm();
    vector<MX> swork(f.getWorkSize());

    // A boolean vector indicated whoch nodes are tainted by substitutions
    vector<bool> tainted(swork.size());

    // Temporary stringstream
    stringstream ss;

    // Construct lookup table for expressions
    std::map<const MXNode*, int> expr_lookup;
    for (int i=0;i<expr.size();++i) {
      expr_lookup[expr[i].operator->()] = i;
    }

    // Construct found map
    std::vector<bool> expr_found(expr.size());

    // Allocate output vector
    vector<MX> f_out(f.getNumOutputs());

    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;

    // expr_lookup iterator
    std::map<const MXNode*, int>::const_iterator it_lookup;

    for (vector<MXAlgEl>::const_iterator it=algorithm.begin(); it!=algorithm.end(); ++it) {

      if (!(it->data).isNull()) {
        // Check if it->data points to a supplied expr
        it_lookup = expr_lookup.find((it->data).operator->());

        if (it->res.front()>=0 && it_lookup!=expr_lookup.end()) {
          // Fill in that expression in-place
          swork[it->res.front()] = exprs[it_lookup->second];
          tainted[it->res.front()] = true;
          expr_found[it_lookup->second] = true;
          continue;
        }
      }

      switch (it->op) {
      case OP_INPUT:
        tainted[it->res.front()] = false;
      case OP_PARAMETER:
        swork[it->res.front()] = it->data;
        tainted[it->res.front()] = false;
        break;
      case OP_OUTPUT:
        f_out[it->res.front()] = swork[it->arg.front()];
        break;
      default:
        {
          bool node_tainted = false;

          input_p.resize(it->arg.size());
          for (int i=0; i<input_p.size(); ++i) {
            int el = it->arg[i];
            if (el>=0) node_tainted =  node_tainted || tainted[el];
            input_p[i] = el<0 ? 0 : &swork[el];
          }

          output_p.resize(it->res.size());
          for (int i=0; i<output_p.size(); ++i) {
            int el = it->res[i];
            output_p[i] = el<0 ? 0 : &swork[el];
            if (el>=0) tainted[el] = node_tainted;
          }

          if (it->res.size()==1 && it->res[0]>=0 && !node_tainted) {
            int el = it->res[0];
            swork[el] = it->data;
          } else {
            const_cast<MX&>(it->data)->evaluateMX(input_p, output_p,
                                                  dummy_p, dummy_p, dummy_p, dummy_p, false);
          }
        }
      }
    }

    bool all_found=true;
    for (int i=0;i<expr.size();++i) {
      all_found = all_found && expr_found[i];
    }

    //casadi_assert_message(all_found,
    //             "MXFunctionInternal::extractNodes(const std::vector<MX>& expr):"
    //             " failed to locate all input expr."
    //             << std::endl << "Here's a boolean list showing which ones where found: "
    //             << expr_found);

    return f_out;

  }

  void extractShared(std::vector<MX>& ex, std::vector<MX>& v, std::vector<MX>& vdef,
                     const std::string& v_prefix, const std::string& v_suffix) {

    // Sort the expression
    MXFunction f(vector<MX>(), ex);
    f.init();

    // Get references to the internal data structures
    const vector<MXAlgEl>& algorithm = f.algorithm();
    vector<MX> work(f.getWorkSize());

    // Count how many times an expression has been used
    vector<int> usecount(work.size(), 0);

    // Remember the origin of every calculation
    vector<pair<int, int> > origin(work.size(), make_pair(-1, -1));

    // Which evaluations to replace
    vector<pair<int, int> > replace;

    // Evaluate the algorithm to identify which evaluations to replace
    int k=0;
    for (vector<MXAlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it, ++k) {
      // Increase usage counters
      switch (it->op) {
      case OP_CONST:
      case OP_PARAMETER:
        break;
      default: // Unary operation, binary operation or output
        for (int c=0; c<it->arg.size(); ++c) {
          if (usecount[it->arg[c]]==0) {
            usecount[it->arg[c]]=1;
          } else if (usecount[it->arg[c]]==1) {
            replace.push_back(origin[it->arg[c]]);
            usecount[it->arg[c]]=-1; // Extracted, do not extract again
          }
        }
      }

      // Perform the operation
      switch (it->op) {
      case OP_OUTPUT:
        break;
      case OP_CONST:
      case OP_PARAMETER:
        usecount[it->res.front()] = -1; // Never extract since it is a primitive type
        break;
      default:
        for (int c=0; c<it->res.size(); ++c) {
          if (it->res[c]>=0) {
            work[it->res[c]] = it->data.getOutput(c);
            usecount[it->res[c]] = 0; // Not (yet) extracted
            origin[it->res[c]] = make_pair(k, c);
          }
        }
        break;
      }
    }

    // New variables and definitions
    v.clear();
    v.reserve(replace.size());
    vdef.clear();
    vdef.reserve(replace.size());

    // Quick return
    if (replace.empty()) return;

    // Sort the elements to be replaced in the order of appearence in the algorithm
    sort(replace.begin(), replace.end());
    vector<pair<int, int> >::const_iterator replace_it=replace.begin();

    // Name of intermediate variables
    stringstream v_name;

    // Arguments for calling the atomic operations
    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;

    // Evaluate the algorithm
    k=0;
    for (vector<MXAlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it, ++k) {
      switch (it->op) {
      case OP_OUTPUT:     ex[it->res.front()] = work[it->arg.front()];      break;
      case OP_CONST:
      case OP_PARAMETER:  work[it->res.front()] = it->data; break;
      default:
        {
          // Pointers to the arguments of the evaluation
          input_p.resize(it->arg.size());
          for (int i=0; i<input_p.size(); ++i) {
            int el = it->arg[i]; // index of the argument
            input_p[i] = el<0 ? 0 : &work[el];
          }

          // Pointers to the result of the evaluation
          output_p.resize(it->res.size());
          for (int i=0; i<output_p.size(); ++i) {
            int el = it->res[i]; // index of the output
            output_p[i] = el<0 ? 0 : &work[el];
          }

          // Evaluate atomic operation
          const_cast<MX&>(it->data)->evaluateMX(input_p, output_p,
                                                dummy_p, dummy_p, dummy_p, dummy_p, false);

          // Possibly replace results with new variables
          for (int c=0; c<it->res.size(); ++c) {
            int ind = it->res[c];
            if (ind>=0 && replace_it->first==k && replace_it->second==c) {
              // Store the result
              vdef.push_back(work[ind]);

              // Create a new variable
              v_name.str(string());
              v_name << v_prefix << v.size() << v_suffix;
              v.push_back(MX::sym(v_name.str()));

              // Use in calculations
              work[ind] = v.back();

              // Go to the next element to be replaced
              replace_it++;
            }
          }
        }
      }
    }
  }

  void printCompact(const MX& ex, std::ostream &stream) {
    // Extract shared subexpressions from ex
    vector<MX> v, vdef;
    vector<MX> ex_extracted(1, ex);
    extractShared(ex_extracted, v, vdef, "@", "");

    // Print the expression without shared subexpressions
    ex_extracted.front().print(stream);

    // Print the shared subexpressions
    if (!v.empty()) {
      stream << endl << "where:" << endl;
      for (int i=0; i<v.size(); ++i) {
        stream << v[i] << " := " << vdef[i] << endl;
      }
    }
  }

  MX jacobian(const MX& ex, const MX &arg) {
    MXFunction temp(arg, ex); // make a runtime
    temp.setOption("name", "helper_jacobian_MX");
    temp.init();
    return temp.jac();
  }

  MX gradient(const MX& ex, const MX &arg) {
    MXFunction temp(arg, ex); // make a runtime
    temp.setOption("name", "helper_gradient_MX");
    temp.init();
    return temp.grad();
  }

  MX tangent(const MX& ex, const MX &arg) {
    MXFunction temp(arg, ex); // make a runtime
    temp.setOption("name", "helper_tangent_MX");
    temp.init();
    return temp.tang();
  }

  MX blockcat(const std::vector< std::vector<MX > > &v) {
    // Quick return if no block rows
    if (v.empty()) return MX::sparse(0, 0);

    // Make sure same number of block columns
    int ncols = v.front().size();
    for (vector<vector<MX> >::const_iterator it=v.begin(); it!=v.end(); ++it) {
      casadi_assert_message(it->size()==ncols, "blockcat: Inconsistent number of blocl columns");
    }

    // Quick return if no block columns
    if (v.front().empty()) return MX::sparse(0, 0);

    // Horizontally concatenate all columns for each row, then vertically concatenate rows
    std::vector<MX> rows;
    for (vector<vector<MX> >::const_iterator it=v.begin(); it!=v.end(); ++it) {
      rows.push_back(horzcat(*it));
    }
    return vertcat(rows);
  }

  MX blockcat(const MX &A, const MX &B, const MX &C, const MX &D) {
    return vertcat(horzcat(A, B), horzcat(C, D));
  }

  MX det(const MX& A) {
    return A->getDeterminant();
  }

  MX inv(const MX& A) {
    return A->getInverse();
  }

  std::vector<MX> getSymbols(const MX& e) {
    MXFunction f(std::vector<MX>(), e);
    f.init();
    return f.getFree();
  }

  std::vector<MX> getSymbols(const std::vector<MX>& e) {
    MXFunction f(std::vector<MX>(), e);
    f.init();
    return f.getFree();
  }

  bool dependsOn(const MX& ex, const std::vector<MX> &arg) {
    if (ex.size()==0) return false;

    // Construct a temporary algorithm
    MXFunction temp(arg, ex);
    temp.init();
    temp.spInit(true);

    for (int i=0;i<temp.getNumInputs();++i) {
      bvec_t* input_ =  get_bvec_t(temp.input(i).data());
      std::fill(input_, input_+temp.input(i).size(), bvec_t(1));
    }
    bvec_t* output_ = get_bvec_t(temp.output().data());
    // Perform a single dependency sweep
    temp.spEvaluate(true);

    // Loop over results
    for (int i=0;i<temp.output().size();++i) {
      if (output_[i]) return true;
    }

    return false;
  }

  MX matrix_expand(const MX& e, const std::vector<MX> &boundary) {
    std::vector<MX> e_v(1, e);
    return matrix_expand(e_v, boundary).at(0);
  }

  std::vector<MX> matrix_expand(const std::vector<MX>& e, const std::vector<MX> &boundary) {

    // Create symbols for boundary nodes
    std::vector<MX> syms(boundary.size());

    for (int i=0;i<syms.size();++i) {
      syms[i] = MX::sym("x", boundary[i].sparsity());
    }

    // Substitute symbols for boundary nodes
    std::vector<MX> ret = graph_substitute(e, boundary, syms);

    // Obtain list of dependents
    std::vector<MX> v = getSymbols(ret);

    // Construct an MXFunction with it
    MXFunction f(v, ret);
    f.init();

    // Expand to SXFunction
    SXFunction s = f.expand();
    s.init();

    return s.call(graph_substitute(v, syms, boundary), true);
  }

  MX kron(const MX& a, const MX& b) {
    const Sparsity &a_sp = a.sparsity();
    MX filler = MX::sparse(b.shape());
    std::vector< std::vector< MX > > blocks(a.size1(), std::vector< MX >(a.size2(), filler));
    for (int i=0;i<a.size1();++i) {
      for (int j=0;j<a.size2();++j) {
        int k = a_sp.getNZ(i, j);
        if (k!=-1) {
          blocks[i][j] = a[k]*b;
        }
      }
    }
    return blockcat(blocks);
  }

  MX solve(const MX& A, const MX& b, const std::string& lsolver, const Dictionary& dict) {
    LinearSolver mysolver(lsolver, A.sparsity(), b.size2());
    mysolver.setOption(dict);
    mysolver.init();
    return mysolver.solve(A, b, false);
  }

  MX pinv(const MX& A, const std::string& lsolver, const Dictionary& dict) {
    if (A.size1()>=A.size2()) {
      return solve(mul(A.T(), A), A.T(), lsolver, dict);
    } else {
      return solve(mul(A, A.T()), A, lsolver, dict).T();
    }
  }

  MX nullspace(const MX& A) {
    SX n = SX::sym("A", A.sparsity());
    SXFunction f(n, nullspace(n));
    f.setOption("name", "nullspace");
    f.init();
    return f(A).at(0);
  }

} // namespace casadi

