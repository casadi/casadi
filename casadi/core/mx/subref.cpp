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

  SubRef* SubRef::clone() const {
    return new SubRef(*this);
  }

  void SubRef::evalD(const cpv_double& input, const pv_double& output,
                         int* itmp, double* rtmp) {
    evalGen<double>(input, output, itmp, rtmp);
  }

  void SubRef::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                          int* itmp, SXElement* rtmp) {
    evalGen<SXElement>(input, output, itmp, rtmp);
  }

  template<typename T>
  void SubRef::evalGen(const std::vector<const T*>& input,
                       const std::vector<T*>& output, int* itmp, T* rtmp) {
    casadi_error("not ready");
  }

  void SubRef::spFwd(const cpv_bvec_t& arg,
                     const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    casadi_error("not ready");
  }

  void SubRef::spAdj(const pv_bvec_t& arg,
                     const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    casadi_error("not ready");
  }

  void SubRef::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "(";
    } else {
      stream << "[" << i_ << ", " << j_ << "])";
    }
  }

  void SubRef::eval(const MXPtrV& input, MXPtrV& output) {
    casadi_error("not ready");
  }

  void SubRef::evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens) {
    casadi_error("not ready");
  }

  void SubRef::evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens) {
    casadi_error("not ready");
  }

  void SubRef::generate(std::ostream &stream, const std::vector<int>& arg,
                                 const std::vector<int>& res, CodeGenerator& gen) const {
    casadi_error("not ready");
  }

} // namespace casadi
