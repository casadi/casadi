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


#include "qpsol.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  // Constructor
  Qpsol::Qpsol(const std::string& name,
               const std::map<std::string, Sparsity> &st)
  : FunctionInternal(name) {
    for (auto i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        A_ = i->second;
      } else if (i->first=="h") {
        H_ = i->second;
      } else {
        casadi_error("Unrecognized field in QP structure: " << i->first);
      }
    }

    n_ = H_.size2();
    nc_ = A_.is_null() ? 0 : A_.size1();

    if (!A_.is_null()) {
      casadi_assert_message(A_.size2()==n_,
        "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :"
        << std::endl <<
        "H: " << H_.dim() << " - A: " << A_.dim() << std::endl <<
        "We need: H_.size2()==A_.size2()" << std::endl);
    }

    casadi_assert_message(H_.is_symmetric(),
      "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
      "H: " << H_.dim() <<
      "We need H square & symmetric" << std::endl);

    // Sparsity
    Sparsity x_sparsity = Sparsity::dense(n_, 1);
    Sparsity bounds_sparsity = Sparsity::dense(nc_, 1);
  }

  Sparsity Qpsol::get_sparsity_in(int ind) const {
    switch (static_cast<QpsolInput>(ind)) {
    case QPSOL_X0:
    case QPSOL_G:
    case QPSOL_LBX:
    case QPSOL_UBX:
    case QPSOL_LAM_X0:
      return get_sparsity_out(QPSOL_X);
    case QPSOL_LBA:
    case QPSOL_UBA:
      return get_sparsity_out(QPSOL_LAM_A);
    case QPSOL_A:
      return A_;
    case QPSOL_H:
      return H_;
    case QPSOL_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Qpsol::get_sparsity_out(int ind) const {
    switch (static_cast<QpsolOutput>(ind)) {
    case QPSOL_COST:
      return Sparsity::scalar();
    case QPSOL_X:
    case QPSOL_LAM_X:
      return Sparsity::dense(n_, 1);
    case QPSOL_LAM_A:
      return Sparsity::dense(nc_, 1);
    case QPSOL_NUM_OUT: break;
    }
    return Sparsity();
  }

  void Qpsol::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);
  }

  Qpsol::~Qpsol() {
  }

  void Qpsol::checkInputs(const double* lbx, const double* ubx,
                          const double* lba, const double* uba) const {
    for (int i=0; i<n_; ++i) {
      double lb = lbx ? lbx[i] : 0., ub = ubx ? ubx[i] : 0.;
      casadi_assert_message(lb <= ub,
                            "LBX[" << i << "] <= UBX[" << i << "] was violated. "
                            << "Got LBX[" << i << "]=" << lb <<
                            " and UBX[" << i << "] = " << ub << ".");
    }
    for (int i=0; i<nc_; ++i) {
      double lb = lba ? lba[i] : 0., ub = uba ? uba[i] : 0.;
      casadi_assert_message(lb <= ub,
                            "LBA[" << i << "] <= UBA[" << i << "] was violated. "
                            << "Got LBA[" << i << "] = " << lb <<
                            " and UBA[" << i << "] = " << ub << ".");
    }
  }

  void Qpsol::generateNativeCode(std::ostream& file) const {
    casadi_error("Qpsol::generateNativeCode not defined for class "
                 << typeid(*this).name());
  }

  std::map<std::string, Qpsol::Plugin> Qpsol::solvers_;

  const std::string Qpsol::infix_ = "qpsol";

  double Qpsol::default_in(int ind) const {
    switch (ind) {
    case QPSOL_LBX:
    case QPSOL_LBA:
      return -std::numeric_limits<double>::infinity();
    case QPSOL_UBX:
    case QPSOL_UBA:
      return std::numeric_limits<double>::infinity();
    default:
      return 0;
    }
  }

} // namespace casadi




