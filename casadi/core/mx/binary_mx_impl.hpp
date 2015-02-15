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


#ifndef CASADI_BINARY_MX_IMPL_HPP
#define CASADI_BINARY_MX_IMPL_HPP

#include "binary_mx.hpp"
#include "mx_tools.hpp"
#include <vector>
#include <sstream>
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../std_vector_tools.hpp"
#include "../casadi_options.hpp"

using namespace std;

namespace casadi {

  template<bool ScX, bool ScY>
  BinaryMX<ScX, ScY>::BinaryMX(Operation op, const MX& x, const MX& y) : op_(op) {
    setDependencies(x, y);
    if (ScX) {
      setSparsity(y.sparsity());
    } else {
      setSparsity(x.sparsity());
    }
  }

  template<bool ScX, bool ScY>
  BinaryMX<ScX, ScY>::~BinaryMX() {
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      casadi_math<double>::printPre(op_, stream);
    } else if (part==1) {
      casadi_math<double>::printSep(op_, stream);
    } else {
      casadi_math<double>::printPost(op_, stream);
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::eval(const cpv_MX& arg, const pv_MX& res) {
    casadi_math<MX>::fun(op_, *arg[0], *arg[1], *res[0]);
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::evalFwd(const std::vector<cpv_MX>& fseed,
                                   const std::vector<pv_MX>& fsens) {
    // Get partial derivatives
    MX pd[2];
    casadi_math<MX>::der(op_, dep(0), dep(1), shared_from_this<MX>(), pd);

    // Propagate forward seeds
    for (int d=0; d<fsens.size(); ++d) {
      *fsens[d][0] = pd[0]*(*fseed[d][0]) + pd[1]*(*fseed[d][1]);
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::evalAdj(const std::vector<pv_MX>& aseed,
                                   const std::vector<pv_MX>& asens) {
    // Get partial derivatives
    MX pd[2];
    casadi_math<MX>::der(op_, dep(0), dep(1), shared_from_this<MX>(), pd);

    // Propagate adjoint seeds
    for (int d=0; d<aseed.size(); ++d) {
      MX s = *aseed[d][0];
      *aseed[d][0] = MX();
      for (int c=0; c<2; ++c) {
        // Get increment of sensitivity c
        MX t = pd[c]*s;

        // If dimension mismatch (i.e. one argument is scalar), then sum all the entries
        if (!t.isScalar() && t.shape() != dep(c).shape()) {
          if (pd[c].shape()!=s.shape()) pd[c] = MX(s.sparsity(), pd[c]);
          t = inner_prod(pd[c], s);
        }

        // Propagate the seeds
        asens[d][c]->addToSum(t);
      }
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::generate(std::ostream &stream,
                                    const std::vector<int>& arg,
                                    const std::vector<int>& res,
                                    CodeGenerator& gen) const {

    // Check if inplace
    bool inplace;
    switch (op_) {
    case OP_ADD:
    case OP_SUB:
    case OP_MUL:
    case OP_DIV:
      inplace = res[0]==arg[0];
    default:
      inplace = false;
    }

    // Print loop and right hand side
    stream << "  for (i=0, "
           << "rr=" << gen.work(res.at(0)) << ", ";
    if (!inplace) stream << "cr=" << gen.work(arg.at(0)) << ", ";
    stream << "cs=" << gen.work(arg.at(1)) << "; i<" << sparsity().nnz() << "; ++i) "
           << "*rr++";

    if (inplace) {
      casadi_math<double>::printSep(op_, stream);
      stream << "=";
      stream << (ScY ? " *cs " : " *cs++ ");
    } else {
      stream << "=";
      casadi_math<double>::printPre(op_, stream);
      stream << (ScX ? " *cr " : " *cr++ ");
      casadi_math<double>::printSep(op_, stream);
      stream << (ScY ? " *cs " : " *cs++ ");
      casadi_math<double>::printPost(op_, stream);
    }

    stream << ";" << endl;
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::evalD(const cpv_double& input, const pv_double& output,
                                     int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                                      int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<bool ScX, bool ScY>
  template<typename T>
  void BinaryMX<ScX, ScY>::evalGen(const std::vector<const T*>& input,
                                   const std::vector<T*>& output, int* itmp, T* rtmp) {
    // Get data
    T* output0 = output[0];
    const T* input0 = input[0];
    const T* input1 = input[1];

    if (!ScX && !ScY) {
      casadi_math<T>::fun(op_, input0, input1, output0, nnz());
    } else if (ScX) {
      casadi_math<T>::fun(op_, *input0, input1, output0, nnz());
    } else {
      casadi_math<T>::fun(op_, input0, *input1, output0, nnz());
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::spFwd(const cpv_bvec_t& arg,
                                 const pv_bvec_t& res,
                                 int* itmp, bvec_t* rtmp) {
    const bvec_t *a0=arg[0], *a1=arg[1];
    bvec_t *r=res[0];
    int n=nnz();
    for (int i=0; i<n; ++i) {
      if (ScX && ScY)
        *r++ = *a0 | *a1;
      else if (ScX && !ScY)
        *r++ = *a0 | *a1++;
      else if (!ScX && ScY)
        *r++ = *a0++ | *a1;
      else
        *r++ = *a0++ | *a1++;
    }
  }

  template<bool ScX, bool ScY>
  void BinaryMX<ScX, ScY>::spAdj(const pv_bvec_t& arg,
                                 const pv_bvec_t& res,
                                 int* itmp, bvec_t* rtmp) {
    bvec_t *a0=arg[0], *a1=arg[1], *r = res[0];
    int n=nnz();
    for (int i=0; i<n; ++i) {
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
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX, ScY>::getUnary(int op) const {
    switch (op_) {
    default: break; // no rule
    }

    // Fallback to default implementation
    return MXNode::getUnary(op);
  }

  template<bool ScX, bool ScY>
  MX BinaryMX<ScX, ScY>::getBinary(int op, const MX& y, bool scX, bool scY) const {
    if (!CasadiOptions::simplification_on_the_fly) return MXNode::getBinary(op, y, scX, scY);

    switch (op_) {
    case OP_ADD:
      if (op==OP_SUB && isEqual(y, dep(0), maxDepth())) return dep(1);
      if (op==OP_SUB && isEqual(y, dep(1), maxDepth())) return dep(0);
      break;
    case OP_SUB:
      if (op==OP_SUB && isEqual(y, dep(0), maxDepth())) return -dep(1);
      if (op==OP_ADD && isEqual(y, dep(1), maxDepth())) return dep(0);
      break;
    default: break; // no rule
    }

    // Fallback to default implementation
    return MXNode::getBinary(op, y, scX, scY);
  }


} // namespace casadi

#endif // CASADI_BINARY_MX_IMPL_HPP
