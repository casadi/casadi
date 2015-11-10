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

    st_.resize(QP_STRUCT_NUM);
    for (std::map<std::string, Sparsity>::const_iterator i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        st_[QP_STRUCT_A]=i->second;
      } else if (i->first=="h") {
        st_[QP_STRUCT_H]=i->second;
      } else {
        casadi_error("Unrecognized field in QP structure: " << i->first);
      }
    }

    addOption("defaults_recipes",    OT_STRINGVECTOR, GenericType(), "",
                                                       "lp", true);

    const Sparsity& A = st_[QP_STRUCT_A];
    const Sparsity& H = st_[QP_STRUCT_H];

    n_ = H.size2();
    nc_ = A.isNull() ? 0 : A.size1();

    if (!A.isNull()) {
      casadi_assert_message(A.size2()==n_,
        "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :"
        << std::endl <<
        "H: " << H.dim() << " - A: " << A.dim() << std::endl <<
        "We need: H.size2()==A.size2()" << std::endl);
    }

    casadi_assert_message(H.is_symmetric(),
      "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
      "H: " << H.dim() <<
      "We need H square & symmetric" << std::endl);

    // Sparsity
    Sparsity x_sparsity = Sparsity::dense(n_, 1);
    Sparsity bounds_sparsity = Sparsity::dense(nc_, 1);

    // Input arguments
    ibuf_.resize(QPSOL_NUM_IN);
    input(QPSOL_X0) = DM::zeros(x_sparsity);
    input(QPSOL_H) = DM::zeros(H);
    input(QPSOL_G) = DM::zeros(x_sparsity);
    input(QPSOL_A) = DM::zeros(A);
    input(QPSOL_LBA) = -DM::inf(bounds_sparsity);
    input(QPSOL_UBA) =  DM::inf(bounds_sparsity);
    input(QPSOL_LBX) = -DM::inf(x_sparsity);
    input(QPSOL_UBX) =  DM::inf(x_sparsity);
    input(QPSOL_LAM_X0) = DM::zeros(x_sparsity);
    //input(QPSOL_LAM_A0) = DM::zeros(x_sparsity);

    // Output arguments
    obuf_.resize(QPSOL_NUM_OUT);
    output(QPSOL_X) = DM::zeros(x_sparsity);
    output(QPSOL_COST) = 0.0;
    output(QPSOL_LAM_X) = DM::zeros(x_sparsity);
    output(QPSOL_LAM_A) = DM::zeros(bounds_sparsity);

    ischeme_ = Function::qpsol_in();
    oscheme_ = Function::qpsol_out();
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
      return st_[QP_STRUCT_A];
    case QPSOL_H:
      return st_[QP_STRUCT_H];
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

  void Qpsol::init() {
    // Call the init method of the base class
    FunctionInternal::init();
  }

  Qpsol::~Qpsol() {
  }

  void Qpsol::checkInputs() const {
    for (int i=0;i<input(QPSOL_LBX).nnz();++i) {
      casadi_assert_message(input(QPSOL_LBX).at(i)<=input(QPSOL_UBX).at(i),
                            "LBX[i] <= UBX[i] was violated for i=" << i
                            << ". Got LBX[i]=" << input(QPSOL_LBX).at(i)
                            << " and UBX[i]=" << input(QPSOL_UBX).at(i));
    }
    for (int i=0;i<input(QPSOL_LBA).nnz();++i) {
      casadi_assert_message(input(QPSOL_LBA).at(i)<=input(QPSOL_UBA).at(i),
                            "LBA[i] <= UBA[i] was violated for i=" << i
                            << ". Got LBA[i]=" << input(QPSOL_LBA).at(i)
                            << " and UBA[i]=" << input(QPSOL_UBA).at(i));
    }
  }

  void Qpsol::generateNativeCode(std::ostream& file) const {
    casadi_error("Qpsol::generateNativeCode not defined for class "
                 << typeid(*this).name());
  }

  std::map<std::string, Qpsol::Plugin> Qpsol::solvers_;

  const std::string Qpsol::infix_ = "qpsol";

  const double& Qpsol::default_in(int ind) const {
    switch (ind) {
    case QPSOL_LBX:
    case QPSOL_LBA:
      return default_minf();
    case QPSOL_UBX:
    case QPSOL_UBA:
      return default_inf();
    default:
      return default_zero();
    }
  }

} // namespace casadi




