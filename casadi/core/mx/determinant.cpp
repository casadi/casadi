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


#include "determinant.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"

using namespace std;

namespace casadi {

  Determinant::Determinant(const MX& x) {
    setDependencies(x);
    setSparsity(Sparsity::dense(1, 1));
  }

  void Determinant::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "det(";
    } else {
      stream << ")";
    }
  }

  void Determinant::eval(const cpv_MX& input, const pv_MX& output) {
    *output[0] = det(*input[0]);
  }

  void Determinant::evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens) {
    const MX& X = dep();
    MX det_X = shared_from_this<MX>();
    MX trans_inv_X = inv(X).T();
    for (int d=0; d<fwdSens.size(); ++d) {
      *fwdSens[d][0] = det_X * inner_prod(trans_inv_X, *fwdSeed[d][0]);
    }
  }

  void Determinant::evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens) {
    const MX& X = dep();
    MX det_X = shared_from_this<MX>();
    MX trans_inv_X = inv(X).T();
    for (int d=0; d<adjSeed.size(); ++d) {
      adjSens[d][0]->addToSum((*adjSeed[d][0]*det_X) * trans_inv_X);
      *adjSeed[d][0] = MX();
    }
  }

} // namespace casadi
