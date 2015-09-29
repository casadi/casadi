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


#include "subref.hpp"

using namespace std;

namespace casadi {

  SubRef::SubRef(const MX& x, const Slice& i, const Slice& j) : i_(i), j_(j) {
    setDependencies(x);
  }

  void SubRef::evalD(const double** arg, double** res, int* iw, double* w) {
    evalGen<double>(arg, res, iw, w);
  }

  void SubRef::evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w) {
    evalGen<SXElement>(arg, res, iw, w);
  }

  template<typename T>
  void SubRef::evalGen(const T* const* arg, T* const* res, int* iw, T* w) {
    casadi_error("not ready");
  }

  void SubRef::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    casadi_error("not ready");
  }

  void SubRef::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    casadi_error("not ready");
  }

  std::string SubRef::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << i_ << ", " << j_ << "]";
    return ss.str();
  }

  void SubRef::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    casadi_error("not ready");
  }

  void SubRef::evalFwd(const std::vector<std::vector<MX> >& fseed,
                       std::vector<std::vector<MX> >& fsens) {
    casadi_error("not ready");
  }

  void SubRef::evalAdj(const std::vector<std::vector<MX> >& aseed,
                       std::vector<std::vector<MX> >& asens) {
    casadi_error("not ready");
  }

  void SubRef::generate(const std::vector<int>& arg, const std::vector<int>& res,
                        CodeGenerator& g) const {
    casadi_error("not ready");
  }

} // namespace casadi
