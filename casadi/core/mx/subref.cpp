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

  void SubRef::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                         std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void SubRef::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                          std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void SubRef::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                           std::vector<T>& rtmp) {
    input[0]->getSub(*output[0], j_, i_);
  }

  void SubRef::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    casadi_error("not ready");
  }

  void SubRef::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "(";
    } else {
      stream << "[" << i_ << ", " << j_ << "])";
    }
  }

  void SubRef::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                          MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                          bool output_given) {
    casadi_error("not ready");
  }

  void SubRef::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                 const std::vector<std::string>& res, CodeGenerator& gen) const {
    casadi_error("not ready");
  }

} // namespace casadi
