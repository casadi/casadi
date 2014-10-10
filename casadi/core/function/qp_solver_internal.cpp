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


#include "qp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
#include <typeinfo>

INPUTSCHEME(QpSolverInput)
OUTPUTSCHEME(QpSolverOutput)

using namespace std;
namespace casadi {

  // Constructor
  QpSolverInternal::QpSolverInternal(const std::vector<Sparsity> &st) : st_(st) {

    casadi_assert_message(st_.size()==QP_STRUCT_NUM, "Problem structure mismatch");

    const Sparsity& A = st_[QP_STRUCT_A];
    const Sparsity& H = st_[QP_STRUCT_H];

    n_ = H.size2();
    nc_ = A.isNull() ? 0 : A.size1();

    if (!A.isNull()) {
      casadi_assert_message(A.size2()==n_,
        "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :"
        << std::endl <<
        "H: " << H.dimString() << " - A: " << A.dimString() << std::endl <<
        "We need: H.size2()==A.size2()" << std::endl);
    }

    casadi_assert_message(H.isSymmetric(),
      "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
      "H: " << H.dimString() <<
      "We need H square & symmetric" << std::endl);

    // Sparsity
    Sparsity x_sparsity = Sparsity::dense(n_, 1);
    Sparsity bounds_sparsity = Sparsity::dense(nc_, 1);

    // Input arguments
    setNumInputs(QP_SOLVER_NUM_IN);
    input(QP_SOLVER_X0) = DMatrix::zeros(x_sparsity);
    input(QP_SOLVER_H) = DMatrix::zeros(H);
    input(QP_SOLVER_G) = DMatrix::zeros(x_sparsity);
    input(QP_SOLVER_A) = DMatrix::zeros(A);
    input(QP_SOLVER_LBA) = -DMatrix::inf(bounds_sparsity);
    input(QP_SOLVER_UBA) =  DMatrix::inf(bounds_sparsity);
    input(QP_SOLVER_LBX) = -DMatrix::inf(x_sparsity);
    input(QP_SOLVER_UBX) =  DMatrix::inf(x_sparsity);
    input(QP_SOLVER_LAM_X0) = DMatrix::zeros(x_sparsity);
    //input(QP_SOLVER_LAM_A0) = DMatrix::zeros(x_sparsity);

    // Output arguments
    setNumOutputs(QP_SOLVER_NUM_OUT);
    output(QP_SOLVER_X) = DMatrix::zeros(x_sparsity);
    output(QP_SOLVER_COST) = 0.0;
    output(QP_SOLVER_LAM_X) = DMatrix::zeros(x_sparsity);
    output(QP_SOLVER_LAM_A) = DMatrix::zeros(bounds_sparsity);

    input_.scheme = SCHEME_QpSolverInput;
    output_.scheme = SCHEME_QpSolverOutput;
  }

  void QpSolverInternal::init() {
    // Call the init method of the base class
    FunctionInternal::init();
  }

  QpSolverInternal::~QpSolverInternal() {
  }

  void QpSolverInternal::evaluate() {
    throw CasadiException("QpSolverInternal::evaluate: Not implemented");
  }

  void QpSolverInternal::solve() {
    throw CasadiException("QpSolverInternal::solve: Not implemented");
  }

  void QpSolverInternal::checkInputs() const {
    for (int i=0;i<input(QP_SOLVER_LBX).size();++i) {
      casadi_assert_message(input(QP_SOLVER_LBX).at(i)<=input(QP_SOLVER_UBX).at(i),
                            "LBX[i] <= UBX[i] was violated for i=" << i
                            << ". Got LBX[i]=" << input(QP_SOLVER_LBX).at(i)
                            << " and UBX[i]=" << input(QP_SOLVER_UBX).at(i));
    }
    for (int i=0;i<input(QP_SOLVER_LBA).size();++i) {
      casadi_assert_message(input(QP_SOLVER_LBA).at(i)<=input(QP_SOLVER_UBA).at(i),
                            "LBA[i] <= UBA[i] was violated for i=" << i
                            << ". Got LBA[i]=" << input(QP_SOLVER_LBA).at(i)
                            << " and UBA[i]=" << input(QP_SOLVER_UBA).at(i));
    }
  }

  void QpSolverInternal::generateNativeCode(std::ostream& file) const {
    casadi_error("QpSolverInternal::generateNativeCode not defined for class "
                 << typeid(*this).name());
  }

  std::map<std::string, QpSolverInternal::Plugin> QpSolverInternal::solvers_;

  const std::string QpSolverInternal::infix_ = "qpsolver";

} // namespace casadi




