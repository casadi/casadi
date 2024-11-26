/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "logsumexp.hpp"
namespace casadi {

  LogSumExp::LogSumExp(const MX& A) {
    set_dep(A);
    set_sparsity(Sparsity::dense(1, 1));
  }

  std::string LogSumExp::disp(const std::vector<std::string>& arg) const {
    return "logsumexp(" + arg.at(0) + ")";
  }

  void LogSumExp::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = logsumexp(arg[0]);
  }

  void LogSumExp::ad_forward(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) const {
    MX max = mmax(dep(0));
    MX expmm = exp(dep(0)-max);
    MX s = sum1(expmm);
    for (casadi_int d=0; d<fsens.size(); ++d) {
      MX v = project(fseed[d][0], dep().sparsity());
      fsens[d][0] = dot(v, expmm)/s;
    }
  }

  void LogSumExp::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) const {
    MX max = mmax(dep(0));
    MX expmm = exp(dep(0)-max);
    MX s = sum1(expmm);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += expmm*aseed[d][0]/s;
    }
  }

  int LogSumExp::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  template<typename T>
  int LogSumExp::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    res[0][0] = casadi_logsumexp(arg[0], dep(0).nnz());
    return 0;
  }

  void LogSumExp::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res,
                        const std::vector<bool>& arg_is_ref,
                        std::vector<bool>& res_is_ref) const {
    // Perform operation inplace
    g << g.workel(res[0]) << " = "
      << g.logsumexp(g.work(arg[0], dep(0).nnz(), arg_is_ref[0]), dep(0).nnz())
      << "\n";
  }

} // namespace casadi
