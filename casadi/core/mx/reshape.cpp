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

using namespace std;

namespace casadi {

  Reshape::Reshape(const MX& x, Sparsity sp) {
    casadi_assert(x.nnz()==sp.nnz());
    setDependencies(x);
    setSparsity(sp);
  }

  Reshape* Reshape::clone() const {
    return new Reshape(*this);
  }

  void Reshape::evaluateD(const double* const* input, double** output,
                          int* itmp, double* rtmp) {
    evaluateGen<double>(input, output, itmp, rtmp);
  }

  void Reshape::evaluateSX(const SXElement* const* input, SXElement** output,
                           int* itmp, SXElement* rtmp) {
    evaluateGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void Reshape::evaluateGen(const T* const* input, T** output, int* itmp, T* rtmp) {
    // Quick return if inplace
    if (input[0]==output[0]) return;

    T* res = output[0];
    const T* arg = input[0];
    copy(arg, arg+nnz(), res);
  }

  void Reshape::propagateSparsity(double** input, double** output, bool fwd) {
    // Quick return if inplace
    if (input[0]==output[0]) return;

    bvec_t *res_ptr = reinterpret_cast<bvec_t*>(output[0]);
    int n = dep().nnz();
    bvec_t *arg_ptr = reinterpret_cast<bvec_t*>(input[0]);
    if (fwd) {
      copy(arg_ptr, arg_ptr+n, res_ptr);
    } else {
      for (int k=0; k<n; ++k) {
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

  void Reshape::eval(const MXPtrV& input, MXPtrV& output) {
    *output[0] = reshape(*input[0], shape());
  }

  void Reshape::evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens) {
    for (int d = 0; d<fwdSens.size(); ++d) {
      *fwdSens[d][0] = reshape(*fwdSeed[d][0], shape());
    }
  }

  void Reshape::evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens) {
    for (int d=0; d<adjSeed.size(); ++d) {
      MX tmp = reshape(*adjSeed[d][0], dep().shape());
      *adjSeed[d][0] = MX();
      adjSens[d][0]->addToSum(tmp);
    }
  }

  void Reshape::generateOperation(std::ostream &stream, const std::vector<int>& arg,
                                  const std::vector<int>& res, CodeGenerator& gen) const {
    if (arg[0]==res[0]) return;
    gen.copyVector(stream, gen.work(arg[0]), nnz(), gen.work(res[0]), "i", false);
  }

  MX Reshape::getReshape(const Sparsity& sp) const {
    return reshape(dep(0), sp);
  }

} // namespace casadi
