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
  CASADI_LINEARSOLVER_LAPACKLU_EXPORT void dgetrf_casadi(int *m, int *n, double *a,
    int *lda, int *ipiv, int *info) {
    dgetrf_(m, n, a, lda, ipiv, info);
  }

  extern "C"
  CASADI_LINEARSOLVER_LAPACKLU_EXPORT void dgetrs_casadi(char* trans, int *n, int *nrhs, double *a,
                          int *lda, int *ipiv, double *b, int *ldb, int *info) {
    dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
  }

  extern "C"
  CASADI_LINEARSOLVER_LAPACKLU_EXPORT void dgeequ_casadi(int *m, int *n, double *a, int *lda,
                          double *r, double *c,
                          double *colcnd, double *rowcnd, double *amax, int *info) {
    dgeequ_(m, n, a, lda, r, c, colcnd, rowcnd, amax, info);
  }

  extern "C"
  CASADI_LINEARSOLVER_LAPACKLU_EXPORT void dlaqge_casadi(int *m, int *n, double *a, int *lda,
                          double *r, double *c,
                          double *colcnd, double *rowcnd, double *amax, char *equed) {
    dlaqge_(m, n, a, lda, r, c, colcnd, rowcnd, amax, equed);
  }

  extern "C"
  int CASADI_LINEARSOLVER_LAPACKLU_EXPORT
  casadi_register_linearsolver_lapacklu(LinearSolverInternal::Plugin* plugin) {
    plugin->creator = LapackLuDense::creator;
    plugin->name = "lapacklu";
    plugin->doc = LapackLuDense::meta_doc.c_str();
    plugin->version = 23;
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
    input(0).get(mat_);

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

        userOut() << "Warning: " << ss.str() << endl;



        if (allow_equilibration_failure_)  userOut() << "Warning: " << ss.str() << endl;
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

  void LapackLuDense::generate(const std::vector<int>& arg, const std::vector<int>& res,
       CodeGenerator& g, int nrhs, bool transpose) const {

    g.declarations << "extern void dlaqge_casadi(int *m, int *n, double *a, int *lda, double *r, ";
    g.declarations << " double *c, double *colcnd, double *rowcnd, double *amax, char *);" << endl;

    g.declarations << "extern void dgetrf_casadi(int *m, int *, double *, int *, int *, int *);";
    g.declarations << endl;

    g.declarations << "extern void dgetrs_casadi(char* trans, int *n, int *nrhs, double *a,";
    g.declarations << "int *lda, int *ipiv, double *b, int *ldb, int *info);" << endl;

    g.declarations << "extern void dgeequ_casadi(int *m, int *n, double *a, int *lda, double *r, "
                      "double *c, double *colcnd, double *rowcnd, double *amax, int *);" << endl;

    g.body << "int ncol = " << ncol() << ";" << endl;
    g.body << "int nrow = " << nrow() << ";" << endl;
    g.body << "double mat[" << ncol()*ncol() <<  "];" << endl;
    g.body << "for (int i=0;i<" << ncol()*nrow() << ";++i) mat[i]=0;" << endl;

    g.body << "double *A = " << g.work(arg[1], input(LINSOL_A).nnz()) << ";" << endl;
    g.body << "const int* colind = " << g.sparsity(input(LINSOL_A).sparsity()) << "+2;" << endl;
    g.body << "const int* row = " << g.sparsity(input(LINSOL_A).sparsity()) << "+" << 3+ncol();
    g.body << ";" << endl;

    g.body << "for (int cc=0; cc<ncol; ++cc) {" << endl;
    g.body << "  for (int el=colind[cc]; el<colind[cc+1]; ++el) {" << endl;
    g.body << "    int rr=row[el];" << endl;
    g.body << "    mat[rr + nrow*cc] = A[el];" << endl;
    g.body << "  }" << endl;
    g.body << "}" << endl;

    g.body << "int ipiv[" << ncol() <<  "];" << endl;

    g.body << "double r[" << ncol() <<  "];" << endl;
    g.body << "double c[" << nrow() <<  "];" << endl;

    g.body << "double colcnd, rowcnd;" << endl;
    g.body << "double amax;" << endl;
    g.body << "int info = -100;" << endl;
    g.body << "dgeequ_casadi(&ncol, &nrow, mat, &ncol, r, ";
    g.body << "c, &colcnd, &rowcnd, &amax, &info);" << endl;
    g.body << "char equed = 'N';" << endl;
    g.body << "if (info!=0) dlaqge_casadi(&ncol, &nrow, mat, &ncol, r, c,";
    g.body << " &colcnd, &rowcnd, &amax, &equed);" << endl;

    g.body << "dgetrf_casadi(&ncol, &ncol, mat, &ncol, ipiv, &info);" << endl;

    g.body << "info = -100;" << endl;
    g.body << "char trans = '" << (transpose ? "T" : "N") << "';" << endl;
    g.body << "int nrhs = " << nrhs << ";" << endl;

    g.body << "double *x = " << g.work(res[0], input(LINSOL_B).nnz()) << ";" << endl;
    g.body << "if (x!=" << g.work(arg[0], input(LINSOL_B).nnz()) << ") {" << endl;
    g.body << "  for (int i=0;i<" <<  input(LINSOL_B).nnz() << ";++i) x[i] = ";
    g.body << g.work(arg[0], input(LINSOL_B).nnz()) << "[i];" << endl;
    g.body << "}" << endl;
    if (transpose) {
      g.body << "if (equed=='C' || equed=='B')" << endl;
      g.body << "  for (int rhs=0; rhs<nrhs; ++rhs)" << endl;
      g.body << "    for (int i=0; i<nrow; ++i)" << endl;
      g.body << "      x[i+rhs*nrow] *= c[i];" << endl;
    } else {
      g.body << "if (equed=='R' || equed=='B')" << endl;
      g.body << "  for (int rhs=0; rhs<nrhs; ++rhs)" << endl;
      g.body << "    for (int i=0; i<ncol; ++i)" << endl;
      g.body << "      x[i+rhs*nrow] *= r[i];" << endl;
    }

    g.body << "    dgetrs_casadi(&trans, &ncol, &nrhs, mat, &ncol, ipiv, x, &ncol, &info);" << endl;

    if (transpose) {
      g.body << "if (equed=='R' || equed=='B')" << endl;
      g.body << "  for (int rhs=0; rhs<nrhs; ++rhs)" << endl;
      g.body << "    for (int i=0; i<ncol; ++i)" << endl;
      g.body << "      x[i+rhs*nrow] *= r[i];" << endl;
    } else {
      g.body << "if (equed=='C' || equed=='B')" << endl;
      g.body << "  for (int rhs=0; rhs<nrhs; ++rhs)" << endl;
      g.body << "    for (int i=0; i<nrow; ++i)" << endl;
      g.body << "      x[i+rhs*nrow] *= c[i];" << endl;
    }

  }

} // namespace casadi
