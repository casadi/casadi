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

  std::string SetSparse::print(const std::vector<std::string>& arg) const {
    if (sparsity().isDense()) {
      return "dense(" + arg.at(0) + ")";
    } else {
      return "set_sparse(" + arg.at(0) + ")";
    }
  }

  template<typename T>
  void SetSparse::evalGen(const T* const* arg, T* const* res,
                          int* iw, T* rtmp) {
    casadi_project(arg[0], dep().sparsity(), res[0], sparsity(), rtmp);
  }

  void SetSparse::evalD(const double** input, double** output,
                            int* iw, double* rtmp) {
    evalGen<double>(input, output, iw, rtmp);
  }

  void SetSparse::evalSX(const SXElement** input, SXElement** output,
                             int* iw, SXElement* rtmp) {
    evalGen<SXElement>(input, output, iw, rtmp);
  }

  void SetSparse::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = arg[0].setSparse(sparsity());
  }

  void SetSparse::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();
    for (int d=0; d<nfwd; ++d) {
      fsens[d][0] = fseed[d][0].setSparse(sparsity(), true);
    }
  }

  void SetSparse::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    int nadj = aseed.size();
    for (int d=0; d<nadj; ++d) {
      asens[d][0] += aseed[d][0].setSparse(dep().sparsity(), true);
    }
  }

  void SetSparse::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* rtmp) {
    sparsity().set(res[0], arg[0], dep().sparsity());
  }

  void SetSparse::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* rtmp) {
    dep().sparsity().bor(arg[0], res[0], sparsity());
    fill(res[0], res[0]+nnz(), 0);
  }

  void SetSparse::generate(const std::vector<int>& arg, const std::vector<int>& res,
                           CodeGenerator& g) const {
    g.body << "  " << g.project(g.work(arg.front(), dep().nnz()), dep(0).sparsity(),
                                g.work(res.front(), nnz()), sparsity(), "w") << endl;
  }


} // namespace casadi

