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

  void InnerProd::eval(const cpv_MX& input, const pv_MX& output) {
    *output[0] = (*input[0])->getInnerProd(*input[1]);
  }

  void InnerProd::evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens) {
    for (int d=0; d<fwdSens.size(); ++d) {
      *fwdSens[d][0] = dep(0)->getInnerProd(*fwdSeed[d][1])
        + (*fwdSeed[d][0])->getInnerProd(dep(1));
    }
  }

  void InnerProd::evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens) {
    for (int d=0; d<adjSeed.size(); ++d) {
      adjSens[d][0]->addToSum(*adjSeed[d][0] * dep(1));
      adjSens[d][1]->addToSum(*adjSeed[d][0] * dep(0));
      *adjSeed[d][0] = MX();
    }
  }

  void InnerProd::evalD(cp_double* input, p_double* output,
                            int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  void InnerProd::evalSX(cp_SXElement* input, p_SXElement* output,
                         int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void InnerProd::evalGen(const T* const* arg, T* const* res,
                          int* itmp, T* rtmp) {
    *res[0] = casadi_dot(dep(0).nnz(), arg[0], 1, arg[1], 1);
  }

  void InnerProd::spFwd(cp_bvec_t* arg,
                        p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    const bvec_t *a0=arg[0], *a1=arg[1];
    bvec_t* r = res[0];
    const int n = dep(0).nnz();
    *r = 0;
    for (int i=0; i<n; ++i) {
      *r |= *a0++ | *a1++;
    }
  }

  void InnerProd::spAdj(p_bvec_t* arg,
                        p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    bvec_t *a0=arg[0], *a1=arg[1], *r=res[0];
    const int n = dep(0).nnz();
    for (int i=0; i<n; ++i) {
      *a0++ |= *r;
      *a1++ |= *r;
    }
    *r = 0;
  }

  void InnerProd::generate(std::ostream &stream, const std::vector<int>& arg,
                                    const std::vector<int>& res, CodeGenerator& gen) const {
    gen.assign(stream, gen.workelement(res[0]),
               gen.casadi_dot(dep().nnz(), gen.work(arg[0]), 1, gen.work(arg[1]), 1));
  }

} // namespace casadi
