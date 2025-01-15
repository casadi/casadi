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


#include "norm.hpp"

namespace casadi {

  Norm::Norm(const MX& x) {
    set_dep(x);
    set_sparsity(Sparsity::scalar());
  }

  std::string NormF::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_F";
  }

  int NormF::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Norm1::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int NormInf::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int NormF::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    eval_gen<SXElem>(arg, res, iw, w);
    return 0;
  }

  int Norm1::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    eval_gen<SXElem>(arg, res, iw, w);
    return 0;
  }

  int NormInf::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    eval_gen<SXElem>(arg, res, iw, w);
    return 0;
  }

  template<typename T>
  int NormF::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    *res[0] = casadi_norm_2(dep().nnz(), arg[0]);
    return 0;
  }

  template<typename T>
  int NormInf::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    *res[0] = casadi_norm_inf(dep().nnz(), arg[0]);
    return 0;
  }

  template<typename T>
  int Norm1::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    *res[0] = casadi_norm_1(dep().nnz(), arg[0]);
    return 0;
  }

  void NormF::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_norm_fro();
  }

  void Norm1::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_norm_1();
  }

  void NormInf::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_norm_inf();
  }

  void Norm2::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_norm_2();
  }

  void NormF::ad_forward(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) const {
    MX self = shared_from_this<MX>();
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = dep(0)->get_dot(fseed[d][0]) / self;
    }
  }

  void NormF::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) const {
    MX self = shared_from_this<MX>();
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += (aseed[d][0]/self) * dep(0);
    }
  }

  void Norm1::ad_forward(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) const {
    MX s = sign(dep(0));
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = s->get_dot(fseed[d][0]);
    }
  }

  void Norm1::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) const {
    MX s = sign(dep(0));
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += s*aseed[d][0];
    }
  }

  void NormInf::ad_forward(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) const {
    MX m = shared_from_this<MX>()==fabs(dep(0));
    MX s = sign(dep(0));
    MX N = sum2(sum1(m));
    for (casadi_int d=0; d<fsens.size(); ++d) {
       fsens[d][0] = dot(s*fseed[d][0], m) / N;
    }
  }

  void NormInf::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) const {
    MX m = shared_from_this<MX>()==fabs(dep(0));
    MX N = sum2(sum1(m));
    MX s = sign(dep(0));
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += (s*aseed[d][0]/N)*m;
    }
  }

  void NormF::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res,
                        const std::vector<bool>& arg_is_ref,
                        std::vector<bool>& res_is_ref) const {
    std::string a = g.work(arg[0], dep(0).nnz(), arg_is_ref[0]);
    g << g.workel(res[0]) << " = " << g.norm_2(dep().nnz(), a) << ";\n";
  }

  void Norm1::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res,
                        const std::vector<bool>& arg_is_ref,
                        std::vector<bool>& res_is_ref) const {
    std::string a = g.work(arg[0], dep(0).nnz(), arg_is_ref[0]);
    g << g.workel(res[0]) << " = " << g.norm_1(dep().nnz(), a) << ";\n";
  }

  void NormInf::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res,
                        const std::vector<bool>& arg_is_ref,
                        std::vector<bool>& res_is_ref) const {
    std::string a = g.work(arg[0], dep(0).nnz(), arg_is_ref[0]);
    g << g.workel(res[0]) << " = " << g.norm_inf(dep().nnz(), a) << ";\n";
  }

  std::string Norm2::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_2";
  }

  std::string Norm1::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_1";
  }

  std::string NormInf::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_inf";
  }

} // namespace casadi
