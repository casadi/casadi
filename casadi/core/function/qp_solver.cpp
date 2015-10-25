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


#include "qp_solver.hpp"
#include <typeinfo>

INPUTSCHEME(QpSolverInput)
OUTPUTSCHEME(QpSolverOutput)

using namespace std;
namespace casadi {

  // Constructor
  QpSolver::QpSolver(const std::string& name,
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

    casadi_assert_message(H.issymmetric(),
      "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
      "H: " << H.dim() <<
      "We need H square & symmetric" << std::endl);

    // Sparsity
    Sparsity x_sparsity = Sparsity::dense(n_, 1);
    Sparsity bounds_sparsity = Sparsity::dense(nc_, 1);

    // Input arguments
    ibuf_.resize(QP_SOLVER_NUM_IN);
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
    obuf_.resize(QP_SOLVER_NUM_OUT);
    output(QP_SOLVER_X) = DMatrix::zeros(x_sparsity);
    output(QP_SOLVER_COST) = 0.0;
    output(QP_SOLVER_LAM_X) = DMatrix::zeros(x_sparsity);
    output(QP_SOLVER_LAM_A) = DMatrix::zeros(bounds_sparsity);

    ischeme_ = IOScheme(SCHEME_QpSolverInput);
    oscheme_ = IOScheme(SCHEME_QpSolverOutput);
  }

  Sparsity QpSolver::get_sparsity_in(int ind) const {
    switch (static_cast<QpSolverInput>(ind)) {
    case QP_SOLVER_X0:
    case QP_SOLVER_G:
    case QP_SOLVER_LBX:
    case QP_SOLVER_UBX:
    case QP_SOLVER_LAM_X0:
      return get_sparsity_out(QP_SOLVER_X);
    case QP_SOLVER_LBA:
    case QP_SOLVER_UBA:
      return get_sparsity_out(QP_SOLVER_LAM_A);
    case QP_SOLVER_A:
      return st_[QP_STRUCT_A];
    case QP_SOLVER_H:
      return st_[QP_STRUCT_H];
    case QP_SOLVER_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity QpSolver::get_sparsity_out(int ind) const {
    switch (static_cast<QpSolverOutput>(ind)) {
    case QP_SOLVER_COST:
      return Sparsity::scalar();
    case QP_SOLVER_X:
    case QP_SOLVER_LAM_X:
      return Sparsity::dense(n_, 1);
    case QP_SOLVER_LAM_A:
      return Sparsity::dense(nc_, 1);
    case QP_SOLVER_NUM_OUT: break;
    }
    return Sparsity();
  }

  void QpSolver::init() {
    // Call the init method of the base class
    FunctionInternal::init();
  }

  QpSolver::~QpSolver() {
  }

  void QpSolver::evaluate() {
    throw CasadiException("QpSolver::evaluate: Not implemented");
  }

  void QpSolver::solve() {
    throw CasadiException("QpSolver::solve: Not implemented");
  }

  void QpSolver::checkInputs() const {
    for (int i=0;i<input(QP_SOLVER_LBX).nnz();++i) {
      casadi_assert_message(input(QP_SOLVER_LBX).at(i)<=input(QP_SOLVER_UBX).at(i),
                            "LBX[i] <= UBX[i] was violated for i=" << i
                            << ". Got LBX[i]=" << input(QP_SOLVER_LBX).at(i)
                            << " and UBX[i]=" << input(QP_SOLVER_UBX).at(i));
    }
    for (int i=0;i<input(QP_SOLVER_LBA).nnz();++i) {
      casadi_assert_message(input(QP_SOLVER_LBA).at(i)<=input(QP_SOLVER_UBA).at(i),
                            "LBA[i] <= UBA[i] was violated for i=" << i
                            << ". Got LBA[i]=" << input(QP_SOLVER_LBA).at(i)
                            << " and UBA[i]=" << input(QP_SOLVER_UBA).at(i));
    }
  }

  void QpSolver::generateNativeCode(std::ostream& file) const {
    casadi_error("QpSolver::generateNativeCode not defined for class "
                 << typeid(*this).name());
  }

  std::map<std::string, QpSolver::Plugin> QpSolver::solvers_;

  const std::string QpSolver::infix_ = "qpsolver";

  const double& QpSolver::default_in(int ind) const {
    switch (ind) {
    case QP_SOLVER_LBX:
    case QP_SOLVER_LBA:
      return default_minf();
    case QP_SOLVER_UBX:
    case QP_SOLVER_UBA:
      return default_inf();
    default:
      return default_zero();
    }
  }

} // namespace casadi




