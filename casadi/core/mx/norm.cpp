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


#include "norm.hpp"
#include "mx_tools.hpp"
#include "../runtime/runtime.hpp"

using namespace std;
namespace casadi {

  Norm::Norm(const MX& x) {
    setDependencies(x);
    setSparsity(Sparsity::scalar());
  }

  void NormF::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_F";
    }
  }

  void NormF::evaluateD(const double* const* input, double** output,
                        int* itmp, double* rtmp) {
    evaluateGen<double>(input, output, itmp, rtmp);
  }

  void NormF::evaluateSX(const SXElement* const* input, SXElement** output,
                         int* itmp, SXElement* rtmp) {
    evaluateGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void NormF::evaluateGen(const T* const* input, T** output, int* itmp, T* rtmp) {
    // Get data
    T* res = output[0];
    const T* arg = input[0];

    // Perform the inner product
    *res = sqrt(casadi_dot(dep().nnz(), arg, 1, arg, 1));
  }

  void NormF::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                         MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                         bool output_given) {
    if (!output_given) {
      *output[0] = (*input[0])->getNormF();
    }

    // Forward sensitivities
    int nfwd = fwdSens.size();
    for (int d=0; d<nfwd; ++d) {
      *fwdSens[d][0] = (*input[0])->getInnerProd(*fwdSeed[d][0]) / (*output[0]);
    }

    // Adjoint sensitivities
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      adjSens[d][0]->addToSum(((*adjSeed[d][0])/(*output[0])) * *input[0]);
      *adjSeed[d][0] = MX();
    }
  }

  void NormF::evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens) {
    MX self = shared_from_this<MX>();
    for (int d=0; d<fwdSens.size(); ++d) {
      *fwdSens[d][0] = dep(0)->getInnerProd(*fwdSeed[d][0]) / self;
    }
  }

  void NormF::evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens) {
    MX self = shared_from_this<MX>();
    for (int d=0; d<adjSeed.size(); ++d) {
      adjSens[d][0]->addToSum(((*adjSeed[d][0])/self) * dep(0));
      *adjSeed[d][0] = MX();
    }
  }

  void NormF::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                const std::vector<std::string>& res, CodeGenerator& gen) const {
    stream << "  *" << res.front() << " = sqrt(" << gen.casadi_dot(dep().nnz(),
                                                                   arg.front(), 1,
                                                                   arg.front(), 1) << ");" << endl;
  }

  void Norm2::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_2";
    }
  }

  void Norm1::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_1";
    }
  }

  void NormInf::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_inf";
    }
  }

} // namespace casadi
