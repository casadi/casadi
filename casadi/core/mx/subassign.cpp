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


#include "subassign.hpp"

using namespace std;

namespace casadi {

  SubAssign::SubAssign(const MX& x, const MX& y, const Slice& i, const Slice& j) : i_(i), j_(j) {
    setDependencies(x, y);
    casadi_error("not ready");
  }

  SubAssign* SubAssign::clone() const {
    return new SubAssign(*this);
  }

  void SubAssign::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output,
                            std::vector<int>& itmp, std::vector<double>& rtmp) {
    evaluateGen<double, DMatrixPtrV, DMatrixPtrVV>(input, output, itmp, rtmp);
  }

  void SubAssign::evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                             std::vector<SXElement>& rtmp) {
    evaluateGen<SXElement, SXPtrV, SXPtrVV>(input, output, itmp, rtmp);
  }

  template<typename T, typename MatV, typename MatVV>
  void SubAssign::evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp,
                              std::vector<T>& rtmp) {
    casadi_error("not ready");
  }

  void SubAssign::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd) {
    casadi_error("not ready");
  }

  void SubAssign::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "(";
    } else if (part==1) {
      stream << "[" << i_ << ", " << j_ << "]=";
    } else {
      stream << ")";
    }
  }

  void SubAssign::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                             MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                             bool output_given) {
    casadi_error("not ready");
  }

  void SubAssign::generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                    const std::vector<std::string>& res, CodeGenerator& gen) const {
    casadi_error("not ready");
  }

} // namespace casadi
