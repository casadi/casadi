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


#include "assertion.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../sx/sx_tools.hpp"

using namespace std;

namespace casadi {

  Assertion::Assertion(const MX& x, const MX& y, const std::string & fail_message)
      : fail_message_(fail_message) {
    casadi_assert_message(y.isScalar(),
                          "Assertion:: assertion expression y must be scalar, but got "
                          << y.dimString());
    setDependencies(x, y);
    setSparsity(x.sparsity());
  }

  void Assertion::printPart(std::ostream &stream, int part) const {
    if (part==0) {
      stream << "assertion(";
    } else if (part==1) {
       stream << ", ";
    } else {
      stream << ")";
    }
  }

  void Assertion::eval(const cpv_MX& input, const pv_MX& output) {
    *output[0] = (*input[0]).attachAssert(*input[1], fail_message_);
  }

  void Assertion::evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens) {
    for (int d=0; d<fwdSens.size(); ++d) {
      *fwdSens[d][0] = *fwdSeed[d][0];
    }
  }

  void Assertion::evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens) {
    for (int d=0; d<adjSeed.size(); ++d) {
      adjSens[d][0]->addToSum(*adjSeed[d][0]);
      *adjSeed[d][0] = MX();
    }
  }

  void Assertion::evalSX(const cpv_SXElement& input, const pv_SXElement& output,
                             int* itmp, SXElement* rtmp) {
    copy(input[0], input[0]+nnz(), output[0]);
  }

  void Assertion::evalD(const cpv_double& input, const pv_double& output,
                            int* itmp, double* rtmp) {
    if (input[1][0]!=1) {
      casadi_error("Assertion error: " << fail_message_);
    }

    copy(input[0], input[0]+nnz(), output[0]);
  }

  void Assertion::spFwd(const cpv_bvec_t& arg,
                        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    copy(arg[0], arg[0]+nnz(), res[0]);
  }

  void Assertion::spAdj(const pv_bvec_t& arg,
                        const pv_bvec_t& res, int* itmp, bvec_t* rtmp) {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    int n = nnz();
    for (int i=0; i<n; ++i) {
      *a++ |= *r;
      *r++ = 0;
    }
  }

} // namespace casadi
