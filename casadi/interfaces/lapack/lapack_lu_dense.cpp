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


#include "lapack_lu_dense.hpp"
#include "../../core/std_vector_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LAPACKLU_EXPORT
  casadi_register_linsol_lapacklu(Linsol::Plugin* plugin) {
    plugin->creator = LapackLuDense::creator;
    plugin->name = "lapacklu";
    plugin->doc = LapackLuDense::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LAPACKLU_EXPORT casadi_load_linsol_lapacklu() {
    Linsol::registerPlugin(casadi_register_linsol_lapacklu);
  }

  LapackLuDense::LapackLuDense(const std::string& name,
                               const Sparsity& sparsity, int nrhs)
    : Linsol(name, sparsity, nrhs) {

    // Equilibrate the matrix
    addOption("equilibration", OT_BOOLEAN, true);
    addOption("allow_equilibration_failure", OT_BOOLEAN, false);
  }

  LapackLuDense::~LapackLuDense() {
  }

  void LapackLuDense::init() {
    // Call the base class initializer
    Linsol::init();

    // Get dimensions
    ncol_ = ncol();
    nrow_ = nrow();

    // Currently only square matrices tested
    if (ncol_!=nrow_) throw CasadiException(
      "LapackLuDense::LapackLuDense: currently only square matrices implemented.");

    // Allocate matrix
    mat_.resize(ncol_*ncol_);
    ipiv_.resize(ncol_);

    // Equilibrate?
    equilibriate_ = option("equilibration").toInt();
    if (equilibriate_) {
      r_.resize(ncol_);
      c_.resize(nrow_);
    }
    equed_ = 'N'; // No equilibration

    // Allow equilibration failures
    allow_equilibration_failure_ = option("allow_equilibration_failure").toInt();
  }

  void LapackLuDense::linsol_factorize(Memory& mem, const double* A) const {
    LapackLuDense& m = const_cast<LapackLuDense&>(*this);

    // Get the elements of the matrix, dense format
    casadi_densify(A, sparsity_, getPtr(m.mat_), false);

    if (equilibriate_) {
      // Calculate the col and row scaling factors
      double colcnd, rowcnd; // ratio of the smallest to the largest col/row scaling factor
      double amax; // absolute value of the largest matrix element
      int info = -100;
      dgeequ_(&m.ncol_, &m.nrow_, getPtr(m.mat_), &m.ncol_, getPtr(m.r_),
              getPtr(m.c_), &colcnd, &rowcnd, &amax, &info);
      if (info < 0)
          throw CasadiException("LapackQrDense::prepare: "
                                "dgeequ_ failed to calculate the scaling factors");
      if (info>0) {
        stringstream ss;
        ss << "LapackLuDense::prepare: ";
        if (info<=ncol_)  ss << (info-1) << "-th row (zero-based) is exactly zero";
        else             ss << (info-1-ncol_) << "-th col (zero-based) is exactly zero";

        userOut() << "Warning: " << ss.str() << endl;



        if (allow_equilibration_failure_)  userOut() << "Warning: " << ss.str() << endl;
        else                              casadi_error(ss.str());
      }

      // Equilibrate the matrix if scaling was successful
      if (info!=0)
        dlaqge_(&m.ncol_, &m.nrow_, getPtr(m.mat_), &m.ncol_, getPtr(m.r_), getPtr(m.c_),
                &colcnd, &rowcnd, &amax, &m.equed_);
      else
        m.equed_ = 'N';
    }

    // Factorize the matrix
    int info = -100;
    dgetrf_(&m.ncol_, &m.ncol_, getPtr(m.mat_), &m.ncol_, getPtr(m.ipiv_), &info);
    casadi_assert_message(info==0, "LapackLuDense::prepare: "
                          "dgetrf_ failed to factorize the Jacobian");
  }

  void LapackLuDense::linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const {
    LapackLuDense& m = const_cast<LapackLuDense&>(*this);

    // Scale the right hand side
    if (tr) {
      m.rowScaling(x, nrhs);
    } else {
      m.colScaling(x, nrhs);
    }

    // Solve the system of equations
    int info = 100;
    char trans = tr ? 'T' : 'N';
    dgetrs_(&trans, &m.ncol_, &nrhs, getPtr(m.mat_), &m.ncol_, getPtr(m.ipiv_), x, &m.ncol_, &info);
    if (info != 0) throw CasadiException("LapackLuDense::solve: "
                                        "failed to solve the linear system");

    // Scale the solution
    if (tr) {
      m.colScaling(x, nrhs);
    } else {
      m.rowScaling(x, nrhs);
    }
  }

  void LapackLuDense::colScaling(double* x, int nrhs) {
    // Scale result if this was done to the matrix
    if (equed_=='R' || equed_=='B')
      for (int rhs=0; rhs<nrhs; ++rhs)
        for (int i=0; i<ncol_; ++i)
          x[i+rhs*nrow_] *= r_[i];
  }

  void LapackLuDense::rowScaling(double* x, int nrhs) {
    // Scale right hand side if this was done to the matrix
    if (equed_=='C' || equed_=='B')
      for (int rhs=0; rhs<nrhs; ++rhs)
        for (int i=0; i<nrow_; ++i)
          x[i+rhs*nrow_] *= c_[i];
  }

} // namespace casadi
