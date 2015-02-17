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


#include "reshape.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/sx_function.hpp"

using namespace std;

namespace casadi {

  Reshape::Reshape(const MX& x, Sparsity sp) {
    casadi_assert(x.nnz()==sp.nnz());
    setDependencies(x);
    setSparsity(sp);
  }

  Reshape* Reshape::clone() const {
    return new Reshape(*this);
  }

  void Reshape::evalD(const cpv_double& input, const pv_double& output,
                          int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  void Reshape::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                           int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void Reshape::evalGen(const std::vector<const T*>& input,
                        const std::vector<T*>& output, int* itmp, T* rtmp) {
    // Quick return if inplace
    if (input[0]==output[0]) return;

    T* res = output[0];
    const T* arg = input[0];
    copy(arg, arg+nnz(), res);
  }

  void Reshape::spFwd(const cpv_bvec_t& arg,
                      const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    copyFwd(arg[0], res[0], nnz());
  }

  void Reshape::spAdj(const pv_bvec_t& arg,
                      const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    copyAdj(arg[0], res[0], nnz());
  }

  void Reshape::printPart(std::ostream &stream, int part) const {
    // For vectors, reshape is also a transpose
    if (dep().isVector(true) && sparsity().isVector(true)) {
      // Print as transpose: X'
      if (part!=0) {
        stream << "'";
      }
    } else {
      // Print as reshape(X) or vec(X)
      if (part==0) {
        if (sparsity().isVector()) {
          stream << "vec(";
        } else {
          stream << "reshape(";
        }
      } else {
        stream << ")";
      }
    }
  }

  void Reshape::eval(const cpv_MX& input, const pv_MX& output) {
    *output[0] = reshape(*input[0], shape());
  }

  void Reshape::evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens) {
    for (int d = 0; d<fwdSens.size(); ++d) {
      *fwdSens[d][0] = reshape(*fwdSeed[d][0], shape());
    }
  }

  void Reshape::evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens) {
    for (int d=0; d<adjSeed.size(); ++d) {
      MX tmp = reshape(*adjSeed[d][0], dep().shape());
      *adjSeed[d][0] = MX();
      adjSens[d][0]->addToSum(tmp);
    }
  }

  void Reshape::generate(std::ostream &stream, const std::vector<int>& arg,
                                  const std::vector<int>& res, CodeGenerator& gen) const {
    if (arg[0]==res[0]) return;
    gen.copyVector(stream, gen.work(arg[0]), nnz(), gen.work(res[0]), "i", false);
  }

  MX Reshape::getReshape(const Sparsity& sp) const {
    return reshape(dep(0), sp);
  }

  MX Reshape::getTranspose() const {
    // For vectors, reshape is also a transpose
    if (dep().isVector(true) && sparsity().isVector(true)) {
      return dep();
    } else {
      return MXNode::getTranspose();
    }
  }

} // namespace casadi
