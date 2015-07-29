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


#include "qcqp_solver_internal.hpp"

INPUTSCHEME(QcqpSolverInput)
OUTPUTSCHEME(QcqpSolverOutput)

using namespace std;
namespace casadi {

  // Constructor
  QcqpSolverInternal::QcqpSolverInternal(const std::map<std::string, Sparsity> &st) {
    st_.resize(QCQP_STRUCT_NUM);
    for (std::map<std::string, Sparsity>::const_iterator i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        st_[QCQP_STRUCT_A]=i->second;
      } else if (i->first=="h") {
        st_[QCQP_STRUCT_H]=i->second;
      } else if (i->first=="p") {
        st_[QCQP_STRUCT_P]=i->second;
      } else {
        casadi_error("Unrecognized field in QCQP structure: " << i->first);
      }
    }

    addOption("defaults_recipes",    OT_STRINGVECTOR, GenericType(), "",
                                                       "qp", true);

    const Sparsity& A = st_[QCQP_STRUCT_A];
    const Sparsity& P = st_[QCQP_STRUCT_P];
    const Sparsity& H = st_[QCQP_STRUCT_H];

    n_ = H.size2();
    nc_ = A.isNull() ? 0 : A.size1();

    if (!A.isNull()) {
      casadi_assert_message(A.size2()==n_,
        "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :"
          << std::endl
          << "H: " << H.dimString() << " - A: " << A.dimString() << std::endl
          << "We need: H.size2()==A.size2()" << std::endl);
    }

    casadi_assert_message(H.size1()==H.size2(),
                          "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
                          "H: " << H.dimString() <<
                          "We need H square & symmetric" << std::endl);

    casadi_assert_message(P.size1()==n_,
                          "Got incompatible dimensions. Number of rows in P (" << P.size1()
                          << ") must match n (" << n_ << ").");

    casadi_assert_message(P.size2() % n_ == 0,
                          "Got incompatible dimensions. Number of cols in P ("
                          << P.size2() << ") must be a multiple of n (" << n_ << ").");


    nq_ = P.size2() / n_;

    // Sparsity
    Sparsity x_sparsity = Sparsity::dense(n_, 1);
    Sparsity bounds_sparsity = Sparsity::dense(nc_, 1);

    // Input arguments
    ibuf_.resize(QCQP_SOLVER_NUM_IN);
    input(QCQP_SOLVER_X0) = DMatrix::zeros(x_sparsity);
    input(QCQP_SOLVER_H) = DMatrix::zeros(H);
    input(QCQP_SOLVER_G) = DMatrix::zeros(x_sparsity);
    input(QCQP_SOLVER_A) = DMatrix::zeros(A);
    input(QCQP_SOLVER_P) = DMatrix::zeros(P);
    input(QCQP_SOLVER_Q) = DMatrix::zeros(nq_*n_, 1);
    input(QCQP_SOLVER_R) = DMatrix::zeros(nq_, 1);
    input(QCQP_SOLVER_LBA) = -DMatrix::inf(bounds_sparsity);
    input(QCQP_SOLVER_UBA) = DMatrix::inf(bounds_sparsity);
    input(QCQP_SOLVER_LBX) = -DMatrix::inf(x_sparsity);
    input(QCQP_SOLVER_UBX) = DMatrix::inf(x_sparsity);

    for (int i=0;i<nq_;++i) {
      DMatrix Pi = input(QCQP_SOLVER_P)(ALL, Slice(i*n_, (i+1)*n_));
      casadi_assert_message(Pi.sparsity().issymmetric(),
                            "We need Pi square & symmetric but got " << Pi.dimString()
                            << " for i = " << i << ".");
    }

    // Output arguments
    obuf_.resize(QCQP_SOLVER_NUM_OUT);
    output(QCQP_SOLVER_X) = DMatrix::zeros(x_sparsity);
    output(QCQP_SOLVER_COST) = 0.0;
    output(QCQP_SOLVER_LAM_X) = DMatrix::zeros(x_sparsity);
    output(QCQP_SOLVER_LAM_A) = DMatrix::zeros(bounds_sparsity);

    ischeme_ = IOScheme(SCHEME_QcqpSolverInput);
    oscheme_ = IOScheme(SCHEME_QcqpSolverOutput);
  }

  void QcqpSolverInternal::init() {
    // Call the init method of the base class
    FunctionInternal::init();
  }

  QcqpSolverInternal::~QcqpSolverInternal() {
  }

  void QcqpSolverInternal::evaluate() {
    throw CasadiException("QcqpSolverInternal::evaluate: Not implemented");
  }

  void QcqpSolverInternal::solve() {
    throw CasadiException("QcqpSolverInternal::solve: Not implemented");
  }

  void QcqpSolverInternal::checkInputs() const {
    for (int i=0;i<input(QCQP_SOLVER_LBX).nnz();++i) {
      casadi_assert_message(input(QCQP_SOLVER_LBX).at(i)<=input(QCQP_SOLVER_UBX).at(i),
                            "LBX[i] <= UBX[i] was violated for i=" << i
                            << ". Got LBX[i]=" << input(QCQP_SOLVER_LBX).at(i)
                            << " and UBX[i]=" << input(QCQP_SOLVER_UBX).at(i));
    }
    for (int i=0;i<input(QCQP_SOLVER_LBA).nnz();++i) {
      casadi_assert_message(input(QCQP_SOLVER_LBA).at(i)<=input(QCQP_SOLVER_UBA).at(i),
                            "LBA[i] <= UBA[i] was violated for i=" << i
                            << ". Got LBA[i]=" << input(QCQP_SOLVER_LBA).at(i)
                            << " and UBA[i]=" << input(QCQP_SOLVER_UBA).at(i));
    }
  }

  std::map<std::string, QcqpSolverInternal::Plugin> QcqpSolverInternal::solvers_;

  const std::string QcqpSolverInternal::infix_ = "qcqpsolver";

} // namespace casadi
