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


#ifndef CASADI_BINARY_MX_IMPL_HPP
#define CASADI_BINARY_MX_IMPL_HPP

#include "binary_mx.hpp"
#include "casadi_misc.hpp"
#include "global_options.hpp"
#include "serializing_stream.hpp"
#include <sstream>
#include <vector>

namespace casadi {

  template<bool ScX, bool ScY>
  BinaryMX<ScX, ScY>::BinaryMX(Operation op, const MX& x, const MX& y) : op_(op) {
    set_dep(x, y);
    if (ScX) {
      set_sparsity(y.sparsity());
    } else {
      set_sparsity(x.sparsity());
    }
  }

  template<bool ScX, bool ScY>
  BinaryMX<ScX, ScY>::~BinaryMX() {
  }

  template<bool ScX, bool ScY>
  std::string BinaryMX<ScX, ScY>::disp(const std::vector<std::string>& arg) const {
    return casadi_math<double>::print(op_, arg.at(0), arg.at(1));
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    casadi_math<MX>::fun(op_, arg[0], arg[1], res[0]);
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                   std::vector<std::vector<MX> >& fsens) const {
    // Get partial derivatives
    MX pd[2];
    casadi_math<MX>::der(op_, dep(0), dep(1), shared_from_this<MX>(), pd);

    // Propagate forward seeds
    for (casadi_int d=0; d<fsens.size(); ++d) {
      if (op_ == OP_IF_ELSE_ZERO) {
        fsens[d][0] = if_else_zero(pd[1], fseed[d][1]);
      } else {
        fsens[d][0] = pd[0]*fseed[d][0] + pd[1]*fseed[d][1];
      }
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                   std::vector<std::vector<MX> >& asens) const {
    // Get partial derivatives
    MX pd[2];
    casadi_math<MX>::der(op_, dep(0), dep(1), shared_from_this<MX>(), pd);

    // Propagate adjoint seeds
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX s = aseed[d][0];
      if (op_ == OP_IF_ELSE_ZERO) {
        // Special case to avoid NaN propagation
        if (!s.is_scalar() && dep(1).is_scalar()) {
          asens[d][1] += dot(dep(0), s);
        } else {
          asens[d][1] += if_else_zero(dep(0), s);
        }
      } else {
        // General case
        for (casadi_int c=0; c<2; ++c) {
          // Get increment of sensitivity c
          MX t = pd[c]*s;

          // If dimension mismatch (i.e. one argument is scalar), then sum all the entries
          if (!t.is_scalar() && t.size() != dep(c).size()) {
            if (pd[c].size()!=s.size()) pd[c] = MX(s.sparsity(), pd[c]);
            t = dot(pd[c], s);
          }

          // Propagate the seeds
          asens[d][c] += t;
        }
      }
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    // Quick return if nothing to do
    if (nnz()==0) return;

    // Check if inplace
    bool inplace;
    switch (op_) {
    case OP_ADD:
    case OP_SUB:
    case OP_MUL:
    case OP_DIV:
      inplace = res[0]==arg[0];
      break;
    default:
      inplace = false;
      break;
    }

    // Scalar names of arguments (start assuming all scalars)
    std::string r = g.workel(res[0]);
    std::string x = g.workel(arg[0]);
    std::string y = g.workel(arg[1]);

    // Avoid emitting '/*' which will be mistaken for a comment
    if (op_==OP_DIV && g.codegen_scalars && dep(1).nnz()==1) {
      y = "(" + y + ")";
    }

    // Codegen loop, if needed
    if (nnz()>1) {
      // Iterate over result
      g.local("rr", "casadi_real", "*");
      g.local("i", "casadi_int");
      g << "for (i=0, " << "rr=" << g.work(res[0], nnz());
      r = "(*rr++)";

      // Iterate over first argument?
      if (!ScX && !inplace) {
        g.local("cr", "const casadi_real", "*");
        g << ", cr=" << g.work(arg[0], dep(0).nnz());
        if (op_==OP_OR || op_==OP_AND) {
          // Avoid short-circuiting with side effects
          x = "cr[i]";
        } else {
          x = "(*cr++)";
        }

      }

      // Iterate over second argument?
      if (!ScY) {
        g.local("cs", "const casadi_real", "*");
        g << ", cs=" << g.work(arg[1], dep(1).nnz());
        if (op_==OP_OR || op_==OP_AND || op_==OP_IF_ELSE_ZERO) {
          // Avoid short-circuiting with side effects
          y = "cs[i]";
        } else {
          y = "(*cs++)";
        }
      }

      // Close loop
      g << "; i<" << nnz() << "; ++i) ";
    }

    // Perform operation
    g << r << " ";
    if (inplace) {
      g << casadi_math<double>::sep(op_) << "= " << y;
    } else {
      g << " = " << g.print_op(op_, x, y);
    }
    g << ";\n";
  }

  template<bool ScX, bool ScY>
  int BinaryMX<ScX, ScY>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  template<bool ScX, bool ScY>
  int BinaryMX<ScX, ScY>::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<bool ScX, bool ScY>
  template<typename T>
  int BinaryMX<ScX, ScY>::
  eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const {
    // Get data
    T* output0 = res[0];
    const T* input0 = arg[0];
    const T* input1 = arg[1];

    if (!ScX && !ScY) {
      casadi_math<T>::fun(op_, input0, input1, output0, nnz());
    } else if (ScX) {
      casadi_math<T>::fun(op_, *input0, input1, output0, nnz());
    } else {
      casadi_math<T>::fun(op_, input0, *input1, output0, nnz());
    }
    return 0;
  }

  template<bool ScX, bool ScY>
  int BinaryMX<ScX, ScY>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a0=arg[0], *a1=arg[1];
    bvec_t *r=res[0];
    casadi_int n=nnz();
    for (casadi_int i=0; i<n; ++i) {
      if (ScX && ScY)
        *r++ = *a0 | *a1;
      else if (ScX && !ScY)
        *r++ = *a0 | *a1++;
      else if (!ScX && ScY)
        *r++ = *a0++ | *a1;
      else
        *r++ = *a0++ | *a1++;
    }
    return 0;
  }

  template<bool ScX, bool ScY>
  int BinaryMX<ScX, ScY>::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a0=arg[0], *a1=arg[1], *r = res[0];
    casadi_int n=nnz();
    for (casadi_int i=0; i<n; ++i) {
      bvec_t s = *r;
      *r++ = 0;
      if (ScX)
        *a0 |= s;
      else
        *a0++ |= s;
      if (ScY)
        *a1 |= s;
      else
        *a1++ |= s;
    }
    return 0;
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX, ScY>::get_unary(casadi_int op) const {
    //switch (op_) {
    //default: break; // no rule
    //}

    // Fallback to default implementation
    return MXNode::get_unary(op);
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX, ScY>::_get_binary(casadi_int op, const MX& y, bool scX, bool scY) const {
    if (!GlobalOptions::simplification_on_the_fly) return MXNode::_get_binary(op, y, scX, scY);

    switch (op_) {
    case OP_ADD:
      if (op==OP_SUB && MX::is_equal(y, dep(0), maxDepth())) return dep(1);
      if (op==OP_SUB && MX::is_equal(y, dep(1), maxDepth())) return dep(0);
      break;
    case OP_SUB:
      if (op==OP_SUB && MX::is_equal(y, dep(0), maxDepth())) return -dep(1);
      if (op==OP_ADD && MX::is_equal(y, dep(1), maxDepth())) return dep(0);
      break;
    default: break; // no rule
    }

    // Fallback to default implementation
    return MXNode::_get_binary(op, y, scX, scY);
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("BinaryMX::op", static_cast<int>(op_));
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    char type_x = ScX;
    char type_y = ScY;
    char type = type_x | (type_y << 1);
    s.pack("BinaryMX::scalar_flags", type);
  }

  template<bool ScX, bool ScY>
  MXNode* BinaryMX<ScX, ScY>::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("BinaryMX::scalar_flags", t);
    bool scX = t & 1;
    bool scY = t & 2;

    if (scX) {
      if (scY) return new BinaryMX<true, true>(s);
      return new BinaryMX<true, false>(s);
    } else {
      if (scY) return new BinaryMX<false, true>(s);
      return new BinaryMX<false, false>(s);
    }
  }

  template<bool ScX, bool ScY>
  BinaryMX<ScX, ScY>::BinaryMX(DeserializingStream& s) : MXNode(s) {
    int op;
    s.unpack("BinaryMX::op", op);
    op_ = Operation(op);
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX, ScY>::get_solve_triu(const MX& r, bool tr) const {
    // Identify systems with the structure (I - R)
    if (!ScX && !ScY && op_ == OP_SUB) {
      // Is the first term a projected unity matrix?
      if (dep(0).is_op(OP_PROJECT) && dep(0).dep(0).is_eye()) {
        // Is the second term strictly lower triangular?
        if (dep(1).is_op(OP_PROJECT) && dep(1).dep(0).sparsity().is_triu(true)) {
          return dep(1).dep(0)->get_solve_triu_unity(r, tr);
        }
      }
    }
    // Fall back to default routine
    return MXNode::get_solve_triu(r, tr);
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX, ScY>::get_solve_tril(const MX& r, bool tr) const {
    // Identify systems with the structure (I - L)
    if (!ScX && !ScY && op_ == OP_SUB) {
      // Is the first term a projected unity matrix?
      if (dep(0).is_op(OP_PROJECT) && dep(0).dep(0).is_eye()) {
        // Is the second term strictly lower triangular?
        if (dep(1).is_op(OP_PROJECT) && dep(1).dep(0).sparsity().is_tril(true)) {
          return dep(1).dep(0)->get_solve_tril_unity(r, tr);
        }
      }
    }
    // Fall back to default routine
    return MXNode::get_solve_tril(r, tr);
  }

} // namespace casadi

#endif // CASADI_BINARY_MX_IMPL_HPP
