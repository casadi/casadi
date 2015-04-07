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


#include "constant_mx.hpp"
#include <vector>
#include <algorithm>
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace casadi {

  ConstantMX::ConstantMX(const Sparsity& sp) {
    setSparsity(sp);
  }

  ConstantMX::~ConstantMX() {
  }

  void ConstantMX::evalD(cp_double* input, p_double* output,
                             int* itmp, double* rtmp) {
  }

  void ConstantMX::evalSX(cp_SXElement* input, p_SXElement* output,
                              int* itmp, SXElement* rtmp) {
  }

  void ConstantMX::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    res[0] = shared_from_this<MX>();
  }

 void ConstantMX::evalFwd(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) {
   MX zero_sens(size1(), size2());
   for (int d=0; d<fsens.size(); ++d) {
     fsens[d][0] = zero_sens;
   }
 }

  void ConstantMX::evalAdj(const std::vector<std::vector<MX> >& aseed,
                           std::vector<std::vector<MX> >& asens) {
  }

  void ConstantMX::spFwd(cp_bvec_t* arg,
                         p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    fill_n(res[0], nnz(), 0);
  }

  void ConstantMX::spAdj(p_bvec_t* arg,
                         p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    fill_n(res[0], nnz(), 0);
  }

  void ConstantDMatrix::generate(std::ostream &stream, const std::vector<int>& arg,
                                          const std::vector<int>& res,
                                          CodeGenerator& gen) const {
    // Print the constant
    int ind = gen.getConstant(x_.data(), true);

    // Copy the constant to the work vector
    gen.copyVector(stream, "c"+gen.numToString(ind), nnz(), gen.work(res.at(0)), "i", false);
  }

  bool ConstantMX::__nonzero__() const {
    if (numel()!=1) casadi_error("Can only determine truth value of scalar MX.");
    if (nnz()!=1) casadi_error("Can only determine truth value of dense scalar MX.");
    return !isZero();
  }

  ConstantMX* ConstantMX::create(const Sparsity& sp, int val) {
    if (sp.isEmpty(true)) {
      return ZeroByZero::getInstance();
    } else {
      switch (val) {
      case 0: return new Constant<CompiletimeConst<0> >(sp);
      case 1: return new Constant<CompiletimeConst<1> >(sp);
      case -1: return new Constant<CompiletimeConst<(-1)> >(sp);
      default: return new Constant<RuntimeConst<int> >(sp, val);
      }
    }
  }

  ConstantMX* ConstantMX::create(const Sparsity& sp, double val) {
    if (sp.isEmpty(true)) {
      return ZeroByZero::getInstance();
    } else {
      int intval(val);
      if (intval-val==0) {
        return create(sp, intval);
      } else {
        return new Constant<RuntimeConst<double> >(sp, val);
      }
    }
  }

  ConstantMX* ConstantMX::create(const Matrix<double>& val) {
    if (val.nnz()==0) {
      return create(val.sparsity(), 0);
    } else if (val.isScalar()) {
      return create(val.sparsity(), val.toScalar());
    } else {
      // Check if all values are the same
      const vector<double> vdata = val.data();
      double v = vdata[0];
      for (vector<double>::const_iterator i=vdata.begin(); i!=vdata.end(); ++i) {
        if (*i!=v) {
          // Values not all the same
          return new ConstantDMatrix(val);
        }
      }

      // All values identical if reached this point
      return create(val.sparsity(), v);
    }
  }

  bool ConstantDMatrix::isZero() const {
    return x_.isZero();
  }

  bool ConstantDMatrix::isOne() const {
    return x_.isOne();
  }

  bool ConstantDMatrix::isMinusOne() const {
    return x_.isMinusOne();
  }

  bool ConstantDMatrix::isIdentity() const {
    return x_.isIdentity();
  }

  // MX ConstantMX::getMultiplication(const MX& y) const {
  //   if (y.isConstant()) {
  //     // Constant folding
  //     DMatrix xv = getMatrixValue();
  //     DMatrix yv = y->getMatrixValue();
  //     return mul(xv, yv);
  //   } else {
  //     return MXNode::getMultiplication(y);
  //   }
  // }

  MX ConstantMX::getInnerProd(const MX& y) const {
    if (y.isConstant()) {
      // Constant folding
      DMatrix xv = getMatrixValue();
      DMatrix yv = y->getMatrixValue();
      return inner_prod(xv, yv);
    } else {
      return MXNode::getInnerProd(y);
    }
  }

  bool ConstantDMatrix::zz_isEqual(const MXNode* node, int depth) const {
    // Check if same node
    const ConstantDMatrix* n = dynamic_cast<const ConstantDMatrix*>(node);
    if (n==0) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (!std::equal(x_.begin(), x_.end(), n->x_.begin())) return false;

    return true;
  }

  void ZeroByZero::printPart(std::ostream &stream, int part) const {
    stream << "0x0";
  }

  MX ZeroByZero::getSetSparse(const Sparsity& sp) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const {
    casadi_assert(nz.empty());
    return MX::zeros(sp);
  }

  MX ZeroByZero::getSetNonzeros(const MX& y, const std::vector<int>& nz) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::getTranspose() const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::getUnary(int op) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::getBinary(int op, const MX& y, bool ScX, bool ScY) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::getReshape(const Sparsity& sp) const {
    casadi_assert(sp.isEmpty());
    return MX::zeros(sp);
  }

} // namespace casadi

