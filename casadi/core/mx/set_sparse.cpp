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


#include "set_sparse.hpp"
#include "mx_tools.hpp"
#include <vector>
#include <sstream>
#include "../std_vector_tools.hpp"

using namespace std;

namespace casadi {

  SetSparse::SetSparse(const MX& x, const Sparsity& sp) {
    setDependencies(x);
    setSparsity(Sparsity(sp));
  }

  SetSparse* SetSparse::clone() const {
    return new SetSparse(*this);
  }

  void SetSparse::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      if (sparsity().isDense()) {
        stream << "dense(";
      } else {
        stream << "set_sparse(";
      }
    } else {
      stream << ")";
    }
  }

  template<typename T>
  void SetSparse::evaluateGen(const T* const* input, T** output, int* itmp, T* rtmp) {
    casadi_project(input[0], dep().sparsity(), output[0], sparsity(), rtmp);
  }

  void SetSparse::evaluateD(const double* const* input, double** output,
                            int* itmp, double* rtmp) {
    evaluateGen<double>(input, output, itmp, rtmp);
  }

  void SetSparse::evaluateSX(const SXElement* const* input, SXElement** output,
                             int* itmp, SXElement* rtmp) {
    evaluateGen<SXElement>(input, output, itmp, rtmp);
  }

  void SetSparse::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                             MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                             bool output_given) {
    // Evaluate function
    if (!output_given) {
      *output[0] = input[0]->setSparse(sparsity());
    }

    // Propagate forward seeds
    int nfwd = fwdSens.size();
    for (int d=0; d<nfwd; ++d) {
      *fwdSens[d][0] = fwdSeed[d][0]->setSparse(sparsity(), true);
    }

    // Propagate adjoint seeds
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      adjSens[d][0]->addToSum(adjSeed[d][0]->setSparse(dep().sparsity(), true));
      *adjSeed[d][0] = MX();
    }
  }

  void SetSparse::evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens) {
    int nfwd = fwdSens.size();
    for (int d=0; d<nfwd; ++d) {
      *fwdSens[d][0] = fwdSeed[d][0]->setSparse(sparsity(), true);
    }
  }

  void SetSparse::evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens) {
    int nadj = adjSeed.size();
    for (int d=0; d<nadj; ++d) {
      adjSens[d][0]->addToSum(adjSeed[d][0]->setSparse(dep().sparsity(), true));
      *adjSeed[d][0] = MX();
    }
  }

  void SetSparse::propagateSparsity(double** input, double** output, bool fwd) {
    bvec_t *inputd = reinterpret_cast<bvec_t*>(input[0]);
    bvec_t *outputd = reinterpret_cast<bvec_t*>(output[0]);
    if (fwd) {
      sparsity().set(outputd, inputd, dep().sparsity());
    } else {
      dep().sparsity().bor(inputd, outputd, sparsity());
      fill(outputd, outputd + nnz(), bvec_t(0));
    }
  }

  void SetSparse::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                    const std::vector<std::string>& res, CodeGenerator& gen) const {
    // Codegen "copy sparse"
    gen.addAuxiliary(CodeGenerator::AUX_PROJECT);

    // Codegen the operation
    int sp_arg = gen.addSparsity(dep(0).sparsity());
    int sp_res = gen.addSparsity(sparsity());
    stream << "  casadi_project(" << arg.front() << ", s" << sp_arg << ", " << res.front()
           << ", s" << sp_res << ", rrr);" << std::endl;
  }


} // namespace casadi

