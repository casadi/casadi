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

#include "../../core/profiling.hpp"
#include "../../core/casadi_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINEARSOLVER_LAPACKLU_EXPORT
  casadi_register_linearsolver_lapacklu(LinearSolverInternal::Plugin* plugin) {
    plugin->creator = LapackLuDense::creator;
    plugin->name = "lapacklu";
    plugin->doc = LapackLuDense::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_LINEARSOLVER_LAPACKLU_EXPORT casadi_load_linearsolver_lapacklu() {
    LinearSolverInternal::registerPlugin(casadi_register_linearsolver_lapacklu);
  }

  LapackLuDense::LapackLuDense(const Sparsity& sparsity, int nrhs)
      : LinearSolverInternal(sparsity, nrhs) {
    // Equilibrate the matrix
    addOption("equilibration", OT_BOOLEAN, true);
    addOption("allow_equilibration_failure", OT_BOOLEAN, false);
  }

  LapackLuDense::~LapackLuDense() {
  }

  void LapackLuDense::init() {
    // Call the base class initializer
    LinearSolverInternal::init();

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
    equilibriate_ = getOption("equilibration").toInt();
    if (equilibriate_) {
      r_.resize(ncol_);
      c_.resize(nrow_);
    }
    equed_ = 'N'; // No equilibration

    // Allow equilibration failures
    allow_equilibration_failure_ = getOption("allow_equilibration_failure").toInt();

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      profileWriteName(CasadiOptions::profilingLog, this, "LapackLUDense",
                       ProfilingData_FunctionType_Other, 2);

      profileWriteSourceLine(CasadiOptions::profilingLog, this, 0, "prepare", -1);
      profileWriteSourceLine(CasadiOptions::profilingLog, this, 1, "solve", -1);
    }
  }

  void LapackLuDense::prepare() {
    double time_start=0;
    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      time_start = getRealTime(); // Start timer
      profileWriteEntry(CasadiOptions::profilingLog, this);
    }
    prepared_ = false;

    // Get the elements of the matrix, dense format
    input(0).get(mat_, DENSE);

    if (equilibriate_) {
      // Calculate the col and row scaling factors
      double colcnd, rowcnd; // ratio of the smallest to the largest col/row scaling factor
      double amax; // absolute value of the largest matrix element
      int info = -100;
      dgeequ_(&ncol_, &nrow_, getPtr(mat_), &ncol_, getPtr(r_),
              getPtr(c_), &colcnd, &rowcnd, &amax, &info);
      if (info < 0)
          throw CasadiException("LapackQrDense::prepare: "
                                "dgeequ_ failed to calculate the scaling factors");
      if (info>0) {
        stringstream ss;
        ss << "LapackLuDense::prepare: ";
        if (info<=ncol_)  ss << (info-1) << "-th row (zero-based) is exactly zero";
        else             ss << (info-1-ncol_) << "-th col (zero-based) is exactly zero";

        cout << "Warning: " << ss.str() << endl;



        if (allow_equilibration_failure_)  cout << "Warning: " << ss.str() << endl;
        else                              casadi_error(ss.str());
      }

      // Equilibrate the matrix if scaling was successful
      if (info!=0)
        dlaqge_(&ncol_, &nrow_, getPtr(mat_), &ncol_, getPtr(r_), getPtr(c_),
                &colcnd, &rowcnd, &amax, &equed_);
      else
        equed_ = 'N';
    }

    // Factorize the matrix
    int info = -100;
    dgetrf_(&ncol_, &ncol_, getPtr(mat_), &ncol_, getPtr(ipiv_), &info);
    if (info != 0) throw CasadiException("LapackLuDense::prepare: "
                                         "dgetrf_ failed to factorize the Jacobian");

    // Success if reached this point
    prepared_ = true;

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      double time_stop = getRealTime(); // Stop timer
      profileWriteTime(CasadiOptions::profilingLog, this, 0, time_stop-time_start,
                       time_stop-time_start);
      profileWriteExit(CasadiOptions::profilingLog, this, time_stop-time_start);
    }
  }

  void LapackLuDense::solve(double* x, int nrhs, bool transpose) {
    double time_start=0;
    if (CasadiOptions::profiling&& CasadiOptions::profilingBinary) {
      time_start = getRealTime(); // Start timer
      profileWriteEntry(CasadiOptions::profilingLog, this);
    }

    // Scale the right hand side
    if (transpose) {
      rowScaling(x, nrhs);
    } else {
      colScaling(x, nrhs);
    }

    // Solve the system of equations
    int info = 100;
    char trans = transpose ? 'T' : 'N';
    dgetrs_(&trans, &ncol_, &nrhs, getPtr(mat_), &ncol_, getPtr(ipiv_), x, &ncol_, &info);
    if (info != 0) throw CasadiException("LapackLuDense::solve: "
                                        "failed to solve the linear system");

    // Scale the solution
    if (transpose) {
      colScaling(x, nrhs);
    } else {
      rowScaling(x, nrhs);
    }

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      double time_stop = getRealTime(); // Stop timer
      profileWriteTime(CasadiOptions::profilingLog, this, 1,
                       time_stop-time_start, time_stop-time_start);
      profileWriteExit(CasadiOptions::profilingLog, this, time_stop-time_start);
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

  LapackLuDense* LapackLuDense::clone() const {
    return new LapackLuDense(*this);
  }

} // namespace casadi
