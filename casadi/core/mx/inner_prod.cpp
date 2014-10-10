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


#include "inner_prod.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../runtime/runtime.hpp"

using namespace std;

namespace casadi {

  InnerProd::InnerProd(const MX& x, const MX& y) {
    casadi_assert(x.sparsity()==y.sparsity());
    setDependencies(x, y);
    setSparsity(Sparsity::scalar());
  }

  void InnerProd::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "inner_prod(";
    } else if (part==1) {
      stream << ", ";
    } else {
      stream << ")";
    }
  }

  void InnerProd::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                             MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                             bool output_given) {
    if (!output_given) {
      *output[0] = (*input[0])->getInnerProd(*input[1]);
    }

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for (int d=0; d<nfwd; ++d) {
      *fwdSens[d][0] = (*input[0])->getInnerProd(*fwdSeed[d][1])
          + (*fwdSeed[d][0])->getInnerProd(*input[1]);
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      adjSens[d][0]->addToSum(*adjSeed[d][0] * *input[1]);
      adjSens[d][1]->addToSum(*adjSeed[d][0] * *input[0]);
      *adjSeed[d][0] = MX();
    }
  }

  void InnerProd::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                            std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void InnerProd::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                             std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void InnerProd::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                              std::vector<T>& rtmp) {
    // Get data
    T& res = output[0]->data().front();
    const vector<T> &arg0 = input[0]->data();
    const vector<T> &arg1 = input[1]->data();
    const int n = arg0.size();

    // Perform the inner product
    res = casadi_dot(n, getPtr(arg0), 1, getPtr(arg1), 1);
  }

  void InnerProd::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    bvec_t& res = *get_bvec_t(output[0]->data());
    bvec_t* arg0 = get_bvec_t(input[0]->data());
    bvec_t* arg1 = get_bvec_t(input[1]->data());
    const int n = input[0]->size();
    if (fwd) {
      res = 0;
      for (int i=0; i<n; ++i) {
        res |= *arg0++ | *arg1++;
      }
    } else {
      for (int i=0; i<n; ++i) {
        *arg0++ |= res;
        *arg1++ |= res;
      }
      res = 0;
    }
  }

  void InnerProd::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                    const std::vector<std::string>& res, CodeGenerator& gen) const {
    stream << "  *" << res.front() << " = "
           << gen.casadi_dot(dep().size(), arg.at(0), 1, arg.at(1), 1) << ";" << endl;
  }

} // namespace casadi
