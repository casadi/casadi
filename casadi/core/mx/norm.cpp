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
#include "mx_tools.hpp"
#include "../runtime/runtime.hpp"

using namespace std;
namespace casadi {

  Norm::Norm(const MX& x) {
    setDependencies(x);
    setSparsity(Sparsity::scalar());
  }

  void NormF::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_F";
    }
  }

  void NormF::evalD(cp_double* input, p_double* output,
                    int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  void NormF::evalSX(cp_SXElement* input, p_SXElement* output,
                         int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void NormF::evalGen(const T* const* arg, T* const* res,
                      int* itmp, T* rtmp) {
    *res[0] = sqrt(casadi_dot(dep().nnz(), arg[0], 1, arg[0], 1));
  }

  void NormF::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = arg[0]->getNormF();
  }

  void NormF::evalFwd(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) {
    MX self = shared_from_this<MX>();
    for (int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = dep(0)->getInnerProd(fseed[d][0]) / self;
    }
  }

  void NormF::evalAdj(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) {
    MX self = shared_from_this<MX>();
    for (int d=0; d<aseed.size(); ++d) {
      asens[d][0] += (aseed[d][0]/self) * dep(0);
    }
  }

  void NormF::generate(std::ostream &stream, const std::vector<int>& arg,
                       const std::vector<int>& res, CodeGenerator& gen) const {
    gen.assign(stream, gen.workelement(res[0]),
               "sqrt(" + gen.casadi_dot(dep().nnz(), gen.work(arg[0]), 1,
                                        gen.work(arg[0]), 1) + ")");
  }

  void Norm2::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_2";
    }
  }

  void Norm1::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_1";
    }
  }

  void NormInf::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "||";
    } else {
      stream << "||_inf";
    }
  }

} // namespace casadi
