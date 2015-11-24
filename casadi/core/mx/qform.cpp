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


#include "qform.hpp"
#include "../runtime/runtime.hpp"

using namespace std;

namespace casadi {

  Qform::Qform(const MX& x, const MX& A) {
    casadi_assert(x.is_dense() && x.is_column());
    casadi_assert(A.is_square() && A.size1()==x.size1());
    setDependencies(x, A);
    setSparsity(Sparsity::scalar());
  }

  std::string Qform::print(const std::vector<std::string>& arg) const {
    return "qform(" + arg.at(0) + ", " + arg.at(1) + ")";
  }

  void Qform::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = qform(arg[0], arg[1]);
  }

  void Qform::evalFwd(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) {
    MX Ax = mul(dep(1), dep(0));
    for (int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = qform(dep(0), fseed[d][1]) + dot(fseed[d][0], Ax);
    }
  }

  void Qform::evalAdj(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) {
    MX Ax = mul(dep(1), dep(0));
    for (int d=0; d<aseed.size(); ++d) {
      asens[d][0] += aseed[d][0] * Ax;
      asens[d][1] += rank1(MX::zeros(dep(1).sparsity()), aseed[d][0], dep(0));
    }
  }

  void Qform::eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    evalGen<double>(arg, res, iw, w, mem);
  }

  void Qform::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem) {
    evalGen<SXElem>(arg, res, iw, w, mem);
  }

  template<typename T>
  void Qform::evalGen(const T** arg, T** res, int* iw, T* w, void* mem) {
    *res[0] = casadi_qform(arg[1], dep(1).sparsity(), arg[0]);
  }

  void Qform::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    /* Result */
    bvec_t &r = *res[0];

    /* Dependency on all x */
    int n = dep(1).size1();
    const bvec_t *x = arg[0];
    for (int i=0; i<n; ++i) r |= *x++;

    /* Dependency on all A */
    int nnz = dep(1).nnz();
    const bvec_t *A = arg[1];
    for (int i=0; i<nnz; ++i) r |= *A++;
  }

  void Qform::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    /* Seed */
    bvec_t &r = *res[0];

    /* Dependency on all x */
    int n = dep(1).size1();
    bvec_t *x = arg[0];
    for (int i=0; i<n; ++i) *x++ |= r;

    /* Dependency on all A */
    int nnz = dep(1).nnz();
    bvec_t *A = arg[1];
    for (int i=0; i<nnz; ++i) *A++ |= r;

    /* Clear seed */
    r = 0;
  }

  void Qform::generate(CodeGenerator& g, const std::string& mem,
                       const std::vector<int>& arg, const std::vector<int>& res) const {
    g.assign(g.body, g.workel(res[0]),
             g.dot(dep().nnz(), g.work(arg[0], dep(0).nnz()),
                   g.work(arg[1], dep(1).nnz())));
  }

} // namespace casadi
