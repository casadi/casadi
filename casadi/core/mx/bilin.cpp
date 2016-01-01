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


#include "bilin.hpp"
#include "../runtime/runtime.hpp"

using namespace std;

namespace casadi {

  Bilin::Bilin(const MX& A, const MX& x, const MX& y) {
    setDependencies(x, A);
    setSparsity(Sparsity::scalar());
  }

  std::string Bilin::print(const std::vector<std::string>& arg) const {
    return "bilin(" + arg.at(0) + ", " + arg.at(1) + ", " + arg.at(2) + ")";
  }

  void Bilin::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = bilin(arg[0], arg[1], arg[2]);
  }

  void Bilin::evalFwd(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) {
    for (int d=0; d<fsens.size(); ++d) {
      fsens[d][0]
        = bilin(fseed[d][0], dep(1), dep(2))
        + bilin(dep(0), fseed[d][1], dep(2))
        + bilin(dep(0), dep(1), fseed[d][2]);
    }
  }

  void Bilin::evalAdj(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) {
    for (int d=0; d<aseed.size(); ++d) {
      asens[d][0] = rank1(project(asens[d][0], sparsity()),
                          aseed[d][0], dep(1), dep(2));
      asens[d][1] += aseed[d][0] * mtimes(dep(0), dep(2));
      asens[d][2] += aseed[d][0] * mtimes(dep(0).T(), dep(1));
    }
  }

  void Bilin::eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    evalGen<double>(arg, res, iw, w, mem);
  }

  void Bilin::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem) {
    evalGen<SXElem>(arg, res, iw, w, mem);
  }

  template<typename T>
  void Bilin::evalGen(const T** arg, T** res, int* iw, T* w, void* mem) {
    *res[0] = casadi_bilin(arg[0], dep(0).sparsity(), arg[1], arg[2]);
  }

  void Bilin::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    /* Get sparsities */
    int ncol_A = sparsity().size2();
    const int *colind_A = dep(0).colind(), *row_A = dep(0).row();

    /* Return value */
    bvec_t r=0;

    /* Loop over the columns of A */
    int cc, rr, el;
    for (cc=0; cc<ncol_A; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
        /* Get the row */
        rr=row_A[el];

        /* Add contribution */
        r |= arg[1][rr] | arg[0][el] | arg[2][cc];
      }
    }
    *res[0] = r;
  }

  void Bilin::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    /* Get sparsities */
    int ncol_A = sparsity().size2();
    const int *colind_A = dep(0).colind(), *row_A = dep(0).row();

    /* Seed */
    bvec_t s=*res[0];
    *res[0] = 0;

    /* Loop over the columns of A */
    int cc, rr, el;
    for (cc=0; cc<ncol_A; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
        /* Get the row */
        rr=row_A[el];

        /* Add contribution */
        arg[0][el] |= s;
        arg[1][rr] |= s;
        arg[2][cc] |= s;
      }
    }
  }

  void Bilin::generate(CodeGenerator& g, const std::string& mem,
                       const std::vector<int>& arg, const std::vector<int>& res) const {
    g.assign(g.body, g.workel(res[0]),
             g.bilin(g.work(arg[0], dep(0).nnz()),
                     sparsity(),
                     g.work(arg[1], dep(1).nnz()),
                     g.work(arg[2], dep(2).nnz())));
  }

} // namespace casadi
