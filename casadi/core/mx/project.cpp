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


#include "project.hpp"
#include <vector>
#include <sstream>
#include "../std_vector_tools.hpp"

using namespace std;

namespace casadi {

  Project::Project(const MX& x, const Sparsity& sp) {
    setDependencies(x);
    setSparsity(Sparsity(sp));
  }

  std::string Project::print(const std::vector<std::string>& arg) const {
    if (sparsity().is_dense()) {
      return "dense(" + arg.at(0) + ")";
    } else {
      return "project(" + arg.at(0) + ")";
    }
  }

  template<typename T>
  void Project::evalGen(const T** arg, T** res, int* iw, T* w, void* mem) {
    casadi_project(arg[0], dep().sparsity(), res[0], sparsity(), w);
  }

  void Project::evalD(const double** arg, double** res, int* iw, double* w, void* mem) {
    evalGen<double>(arg, res, iw, w, mem);
  }

  void Project::evalSX(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem) {
    evalGen<SXElem>(arg, res, iw, w, mem);
  }

  void Project::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = project(arg[0], sparsity());
  }

  void Project::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
    int nfwd = fsens.size();
    for (int d=0; d<nfwd; ++d) {
      fsens[d][0] = project(fseed[d][0], sparsity(), true);
    }
  }

  void Project::evalAdj(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) {
    int nadj = aseed.size();
    for (int d=0; d<nadj; ++d) {
      asens[d][0] += project(aseed[d][0], dep().sparsity(), true);
    }
  }

  void Project::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    sparsity().set(res[0], arg[0], dep().sparsity());
  }

  void Project::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    dep().sparsity().bor(arg[0], res[0], sparsity());
    fill(res[0], res[0]+nnz(), 0);
  }

  void Project::generate(CodeGenerator& g, const std::string& mem,
                         const std::vector<int>& arg, const std::vector<int>& res) const {
    g.body << "  " << g.project(g.work(arg.front(), dep().nnz()), dep(0).sparsity(),
                                g.work(res.front(), nnz()), sparsity(), "w") << endl;
  }


} // namespace casadi

