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


#include "lapack_lu.hpp"
#include "../../core/std_vector_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LAPACKLU_EXPORT
  casadi_register_linsol_lapacklu(Linsol::Plugin* plugin) {
    plugin->creator = LapackLu::creator;
    plugin->name = "lapacklu";
    plugin->doc = LapackLu::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LAPACKLU_EXPORT casadi_load_linsol_lapacklu() {
    Linsol::registerPlugin(casadi_register_linsol_lapacklu);
  }

  LapackLu::LapackLu(const std::string& name,
                               const Sparsity& sparsity, int nrhs)
    : Linsol(name, sparsity, nrhs) {

    // Equilibrate the matrix
    addOption("equilibration", OT_BOOL, true);
    addOption("allow_equilibration_failure", OT_BOOL, false);
  }

  LapackLu::~LapackLu() {
  }

  void LapackLu::init() {
    // Call the base class initializer
    Linsol::init();

    // Get dimensions
    casadi_assert(nrow()==ncol());

    // Equilibrate?
    equilibriate_ = option("equilibration").to_int();

    // Allow equilibration failures
    allow_equilibration_failure_ = option("allow_equilibration_failure").to_int();
  }

  Memory* LapackLu::memory() const {
    LapackLuMemory* m = new LapackLuMemory();
    try {
      // Allocate matrix
      m->mat.resize(nrow() * ncol());
      m->ipiv.resize(ncol());

      // Equilibration
      if (equilibriate_) {
        m->r.resize(nrow());
        m->c.resize(ncol());
      }
      m->equed = 'N'; // No equilibration

      return m;
    } catch (...) {
      delete m;
      return 0;
    }
  }

  void LapackLu::linsol_factorize(Memory& mem, const double* A) const {
    LapackLuMemory& m = dynamic_cast<LapackLuMemory&>(mem);

    // Dimensions
    int nrow = this->nrow();
    int ncol = this->ncol();

    // Get the elements of the matrix, dense format
    casadi_densify(A, sparsity_, get_ptr(m.mat), false);

    if (equilibriate_) {
      // Calculate the col and row scaling factors
      double colcnd, rowcnd; // ratio of the smallest to the largest col/row scaling factor
      double amax; // absolute value of the largest matrix element
      int info = -100;
      dgeequ_(&ncol, &nrow, get_ptr(m.mat), &ncol, get_ptr(m.r),
              get_ptr(m.c), &colcnd, &rowcnd, &amax, &info);
      if (info < 0)
          throw CasadiException("LapackQrDense::prepare: "
                                "dgeequ_ failed to calculate the scaling factors");
      if (info>0) {
        stringstream ss;
        ss << "LapackLu::prepare: ";
        if (info<=ncol)  ss << (info-1) << "-th row (zero-based) is exactly zero";
        else             ss << (info-1-ncol) << "-th col (zero-based) is exactly zero";

        userOut() << "Warning: " << ss.str() << endl;



        if (allow_equilibration_failure_)  userOut() << "Warning: " << ss.str() << endl;
        else                              casadi_error(ss.str());
      }

      // Equilibrate the matrix if scaling was successful
      if (info!=0)
        dlaqge_(&ncol, &nrow, get_ptr(m.mat), &ncol, get_ptr(m.r), get_ptr(m.c),
                &colcnd, &rowcnd, &amax, &m.equed);
      else
        m.equed = 'N';
    }

    // Factorize the matrix
    int info = -100;
    dgetrf_(&ncol, &ncol, get_ptr(m.mat), &ncol, get_ptr(m.ipiv), &info);
    casadi_assert_message(info==0, "LapackLu::prepare: "
                          "dgetrf_ failed to factorize the Jacobian");
  }

  void LapackLu::linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const {
    LapackLuMemory& m = dynamic_cast<LapackLuMemory&>(mem);

    // Dimensions
    int nrow = this->nrow();
    int ncol = this->ncol();

    // Scale the right hand side
    if (tr) {
      if (m.equed=='C' || m.equed=='B')
        for (int rhs=0; rhs<nrhs; ++rhs)
          for (int i=0; i<nrow; ++i)
            x[i+rhs*nrow] *= m.c[i];
    } else {
      if (m.equed=='R' || m.equed=='B')
        for (int rhs=0; rhs<nrhs; ++rhs)
          for (int i=0; i<ncol; ++i)
            x[i+rhs*nrow] *= m.r[i];
    }

    // Solve the system of equations
    int info = 100;
    char trans = tr ? 'T' : 'N';
    dgetrs_(&trans, &ncol, &nrhs, get_ptr(m.mat), &ncol, get_ptr(m.ipiv), x, &ncol, &info);
    if (info != 0) throw CasadiException("LapackLu::solve: "
                                        "failed to solve the linear system");

    // Scale the solution
    if (tr) {
      if (m.equed=='R' || m.equed=='B')
        for (int rhs=0; rhs<nrhs; ++rhs)
          for (int i=0; i<ncol; ++i)
            x[i+rhs*nrow] *= m.r[i];
    } else {
      if (m.equed=='C' || m.equed=='B')
        for (int rhs=0; rhs<nrhs; ++rhs)
          for (int i=0; i<nrow; ++i)
            x[i+rhs*nrow] *= m.c[i];
    }
  }

} // namespace casadi
