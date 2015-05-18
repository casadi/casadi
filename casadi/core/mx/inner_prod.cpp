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

  std::string InnerProd::print(const std::vector<std::string>& arg) const {
    return "inner_prod(" + arg.at(0) + ", " + arg.at(1) + ")";
  }

  void InnerProd::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = arg[0]->getInnerProd(arg[1]);
  }

  void InnerProd::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    for (int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = dep(0)->getInnerProd(fseed[d][1])
        + fseed[d][0]->getInnerProd(dep(1));
    }
  }

  void InnerProd::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    for (int d=0; d<aseed.size(); ++d) {
      asens[d][0] += aseed[d][0] * dep(1);
      asens[d][1] += aseed[d][0] * dep(0);
    }
  }

  void InnerProd::evalD(const double** input, double** output,
                            int* iw, double* w) {
    evalGen<double>(input, output, iw, w);
  }

  void InnerProd::evalSX(const SXElement** input, SXElement** output,
                         int* iw, SXElement* w) {
    evalGen<SXElement>(input, output, iw, w);
  }

  template<typename T>
  void InnerProd::evalGen(const T* const* arg, T* const* res,
                          int* iw, T* w) {
    *res[0] = casadi_inner_prod(dep(0).nnz(), arg[0], arg[1]);
  }

  void InnerProd::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    const bvec_t *a0=arg[0], *a1=arg[1];
    bvec_t* r = res[0];
    const int n = dep(0).nnz();
    *r = 0;
    for (int i=0; i<n; ++i) {
      *r |= *a0++ | *a1++;
    }
  }

  void InnerProd::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    bvec_t *a0=arg[0], *a1=arg[1], *r=res[0];
    const int n = dep(0).nnz();
    for (int i=0; i<n; ++i) {
      *a0++ |= *r;
      *a1++ |= *r;
    }
    *r = 0;
  }

  void InnerProd::generate(const std::vector<int>& arg, const std::vector<int>& res,
                           CodeGenerator& g) const {
    g.assign(g.body, g.workel(res[0], nnz()),
               g.inner_prod(dep().nnz(), g.work(arg[0], dep(0).nnz()),
                              g.work(arg[1], dep(1).nnz())));
  }

} // namespace casadi
