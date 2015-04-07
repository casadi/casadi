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


#include "inverse.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"

using namespace std;

namespace casadi {

  Inverse::Inverse(const MX& x) {
    casadi_assert_message(x.size1()==x.size2(),
                          "Inverse: matrix must be square, but you supllied " << x.dimString());
    setDependencies(x);
    setSparsity(Sparsity::dense(x.size1(), x.size2()));
  }

  void Inverse::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "inv(";
    } else {
      stream << ")";
    }
  }

  void Inverse::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = inv(arg[0]);
  }

  void Inverse::evalFwd(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) {
    MX inv_X = shared_from_this<MX>();
    for (int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = -mul(inv_X, mul(fseed[d][0], inv_X));
    }
  }

  void Inverse::evalAdj(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) {
    MX inv_X = shared_from_this<MX>();
    MX trans_inv_X = inv_X.T();
    for (int d=0; d<aseed.size(); ++d) {
      asens[d][0] -= mul(trans_inv_X, mul(aseed[d][0], trans_inv_X));
    }
  }

} // namespace casadi
