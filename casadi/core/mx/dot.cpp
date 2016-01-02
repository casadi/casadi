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


#include "dot.hpp"
#include "../runtime/runtime.hpp"

using namespace std;

namespace casadi {

  Dot::Dot(const MX& x, const MX& y) {
    casadi_assert(x.sparsity()==y.sparsity());
    setDependencies(x, y);
    setSparsity(Sparsity::scalar());
  }

  std::string Dot::print(const std::vector<std::string>& arg) const {
    return "dot(" + arg.at(0) + ", " + arg.at(1) + ")";
  }

  void Dot::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = arg[0]->getDot(arg[1]);
  }

  void Dot::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    for (int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = dep(0)->getDot(fseed[d][1])
        + fseed[d][0]->getDot(dep(1));
    }
  }

  void Dot::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    for (int d=0; d<aseed.size(); ++d) {
      asens[d][0] += aseed[d][0] * dep(1);
      asens[d][1] += aseed[d][0] * dep(0);
    }
  }

  void Dot::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    evalGen<double>(arg, res, iw, w, mem);
  }

  void Dot::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w, mem);
  }

  template<typename T>
  void Dot::evalGen(const T** arg, T** res, int* iw, T* w, int mem) const {
    *res[0] = casadi_dot(dep(0).nnz(), arg[0], arg[1]);
  }

  void Dot::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    const bvec_t *a0=arg[0], *a1=arg[1];
    bvec_t* r = res[0];
    const int n = dep(0).nnz();
    *r = 0;
    for (int i=0; i<n; ++i) {
      *r |= *a0++ | *a1++;
    }
  }

  void Dot::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    bvec_t *a0=arg[0], *a1=arg[1], *r=res[0];
    const int n = dep(0).nnz();
    for (int i=0; i<n; ++i) {
      *a0++ |= *r;
      *a1++ |= *r;
    }
    *r = 0;
  }

  void Dot::generate(CodeGenerator& g, const std::string& mem,
                           const std::vector<int>& arg, const std::vector<int>& res) const {
    g.assign(g.body, g.workel(res[0]),
             g.dot(dep().nnz(), g.work(arg[0], dep(0).nnz()),
                   g.work(arg[1], dep(1).nnz())));
  }

} // namespace casadi
