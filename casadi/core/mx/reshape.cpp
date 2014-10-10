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


#include "reshape.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/sx_function.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace casadi {

  Reshape::Reshape(const MX& x, Sparsity sp) {
    casadi_assert(x.size()==sp.size());
    setDependencies(x);
    setSparsity(sp);
  }

  Reshape* Reshape::clone() const {
    return new Reshape(*this);
  }

  void Reshape::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                          std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void Reshape::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                           std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void Reshape::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                            std::vector<T>& rtmp) {
    // Quick return if inplace
    if (input[0]==output[0]) return;

    vector<T>& res = output[0]->data();
    const vector<T>& arg = input[0]->data();
    copy(arg.begin(), arg.end(), res.begin());
  }

  void Reshape::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    // Quick return if inplace
    if (input[0]==output[0]) return;

    bvec_t *res_ptr = get_bvec_t(output[0]->data());
    vector<double>& arg = input[0]->data();
    bvec_t *arg_ptr = get_bvec_t(arg);
    if (fwd) {
      copy(arg_ptr, arg_ptr+arg.size(), res_ptr);
    } else {
      for (int k=0; k<arg.size(); ++k) {
        *arg_ptr++ |= *res_ptr;
        *res_ptr++ = 0;
      }
    }
  }

  void Reshape::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "reshape(";
    } else {
      stream << ")";
    }
  }

  void Reshape::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                           MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                           bool output_given) {
    // Quick return if inplace
    if (input[0]==output[0]) return;

    if (!output_given) {
      *output[0] = reshape(*input[0], shape());
    }

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for (int d = 0; d<nfwd; ++d) {
      *fwdSens[d][0] = reshape(*fwdSeed[d][0], shape());
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      MX& aseed = *adjSeed[d][0];
      MX& asens = *adjSens[d][0];
      asens.addToSum(reshape(aseed, dep().shape()));
      aseed = MX();
    }
  }

  void Reshape::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                  const std::vector<std::string>& res, CodeGenerator& gen) const {
    // Quick return if inplace
    if (arg[0].compare(res[0])==0) return;

    stream << "  for (i=0; i<" << size() << "; ++i) " << res.front()
           << "[i] = " << arg.front() << "[i];" << endl;
  }

  MX Reshape::getReshape(const Sparsity& sp) const {
    return reshape(dep(0), sp);
  }

} // namespace casadi
