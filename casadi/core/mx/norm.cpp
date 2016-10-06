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
#include "../runtime/runtime.hpp"

using namespace std;
namespace casadi {

  Norm::Norm(const MX& x) {
    setDependencies(x);
    setSparsity(Sparsity::scalar());
  }

  std::string NormF::print(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_F";
  }

  void NormF::eval(const double** arg, double** res, int* iw, double* w, int mem) const {
    evalGen<double>(arg, res, iw, w, mem);
  }

  void NormF::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w, mem);
  }

  template<typename T>
  void NormF::evalGen(const T** arg, T** res, int* iw, T* w, int mem) const {
    *res[0] = casadi_norm_2(dep().nnz(), arg[0]);
  }

  void NormF::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = arg[0]->getNormF();
  }

  void NormF::evalFwd(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) {
    MX self = shared_from_this<MX>();
    for (int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = dep(0)->getDot(fseed[d][0]) / self;
    }
  }

  void NormF::evalAdj(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) {
    MX self = shared_from_this<MX>();
    for (int d=0; d<aseed.size(); ++d) {
      asens[d][0] += (aseed[d][0]/self) * dep(0);
    }
  }

  void NormF::generate(CodeGenerator& g, const std::string& mem,
                       const std::vector<int>& arg, const std::vector<int>& res) const {
    g.assign(g.body, g.workel(res[0]),
             "sqrt(" + g.dot(dep().nnz(), g.work(arg[0], dep(0).nnz()),
                             g.work(arg[0], dep(0).nnz())) + ")");
  }

  std::string Norm2::print(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_2";
  }

  std::string Norm1::print(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_1";
  }

  std::string NormInf::print(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_inf";
  }

} // namespace casadi
