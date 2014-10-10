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


#include "slicot_tools.hpp"
#include <cassert>

/// \cond INTERNAL

// Need an 8-byte integer since libslicot0 is compiled with  fdefault-integer-8
typedef long long int f_int;

extern "C" {
  int mb03vd_(f_int* n, f_int* p, f_int* ilo, f_int* ihi, double *a, f_int* lda1, f_int* lda2,
              double* tau, f_int* ldtau, double* dwork, f_int *info);
  int mb03vy_(f_int* n, f_int* p, f_int* ilo, f_int* ihi, double *a, f_int* lda1, f_int* lda2,
              const double* tau, f_int* ldtau, double* dwork, f_int *ld_work, f_int *info);
  int mb03wd_(char* job, char* compz, f_int* n, f_int* p, f_int* ilo, f_int* ihi, f_int* iloz,
              f_int* ihiz, double *h, f_int* ldh1, f_int* ldh2, double* z, f_int* ldz1,
              f_int* ldz2, double* wr, double *wi, double* dwork, f_int *ld_work, f_int *info);
}

namespace casadi {
  void slicot_mb03vd(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2, double * tau,
                     int ldtau, double * dwork) {
     if (dwork==0) {
       std::vector<double> work = std::vector<double>(n);
       slicot_mb03vd(n, p, ilo, ihi, a, lda1, lda2, tau, ldtau, &work[0]);
       return;
     }
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int lda1_ = lda1;
     f_int lda2_ = lda2;
     f_int ldtau_ = ldtau;
     f_int ret_ = 0;

     mb03vd_(&n_, &p_, &ilo_, &ihi_, a, &lda1_, &lda2_, tau, &ldtau_, dwork, &ret_);

     if (ret_<0) {
       casadi_error("mb03vd wrong arguments:" << ret_);
     } else if (ret_>0) {
       casadi_error("mb03vd error code:" << ret_);
     }


  }

  void slicot_mb03vy(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2,
                     const double * tau, int ldtau, double * dwork, int ldwork) {
     if (dwork==0) {
       std::vector<double> work = std::vector<double>(4*n);
       slicot_mb03vy(n, p, ilo, ihi, a, lda1, lda2, tau, ldtau, &work[0], 4*n);
       return;
     }
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int lda1_ = lda1;
     f_int lda2_ = lda2;
     f_int ldtau_ = ldtau;
     f_int ldwork_ = ldwork;
     f_int ret_=0;
     mb03vy_(&n_, &p_, &ilo_, &ihi_, a, &lda1_, &lda2_, tau, &ldtau_, dwork, &ldwork_, &ret_);

     if (ret_ < 0) {
       casadi_error("mb03vy wrong arguments:" << ret_);
     } else if (ret_>0) {
       casadi_error("mb03vy error code:" << ret_);
     }


  }

  void slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi, int iloz, int ihiz,
                     double *h, int ldh1, int ldh2, double* z, int ldz1, int ldz2, double* wr,
                     double *wi, double * dwork, int ldwork) {
if (dwork==0) {
       std::vector<double> work = std::vector<double>(ihi-ilo+p-1);
       slicot_mb03wd(job, compz, n, p, ilo, ihi, iloz, ihiz, h, ldh1, ldh2, z, ldz1, ldz2, wr, wi,
                     &work[0], ihi-ilo+p-1);
       return;
     }
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int iloz_ = ilo;
     f_int ihiz_ = ihi;
     f_int ldh1_ = ldh1;
     f_int ldh2_ = ldh2;
     f_int ldz1_ = ldz1;
     f_int ldz2_ = ldz2;
     f_int ldwork_ = ldwork;
     f_int ret_ = 0;
     mb03wd_(&job, &compz, &n_, &p_, &ilo_, &ihi_, &iloz_, &ihiz_, h, &ldh1_, &ldh2_,
             z, &ldz1_, &ldz2_, wr, wi, dwork, &ldwork_, &ret_);

     if (ret_<0) {
       casadi_error("mb03wd wrong arguments:" << ret_);
     } else if (ret_>0) {
       casadi_error("mb03wd error code:" << ret_);
     }
  }

   void slicot_periodic_schur(int n, int K, const std::vector< double > & a,
                              std::vector< double > & t,  std::vector< double > & z,
                              std::vector<double> &eig_real, std::vector<double> &eig_imag,
                              double num_zero) {
     std::vector<double> dwork(std::max(n+K-2, 4*n)+(n-1)*K);
     slicot_periodic_schur(n, K, a, t, z, dwork, eig_real, eig_imag, num_zero);
   }

   void slicot_periodic_schur(int n, int K, const std::vector< double > & a,
                              std::vector< double > & t,  std::vector< double > & z,
                              std::vector<double> &dwork, std::vector<double> &eig_real,
                              std::vector<double> &eig_imag, double num_zero) {
    int mem_base = std::max(n+K-2, 4*n);
    int mem_needed = mem_base+(n-1)*K;

    if (eig_real.size()!=n) {
      eig_real.resize(n);
    }

    if (eig_imag.size()!=n) {
      eig_imag.resize(n);
    }

    if (dwork.size()==0) {
      dwork.resize(mem_needed);
    } else if (dwork.size()<mem_needed) {
      casadi_warning("You provided working memory for slicot_periodic_schur, but it is "
                     "insufficient. We need size " << mem_needed
                     << " but supplied size is " << dwork.size() <<".");
      dwork.resize(mem_needed);
    } else {
      mem_needed = dwork.size();
    }

    t.resize(n*n*K);

    // a is immutable, we need a mutable pointer, so we use available buffer
    z = a;

    slicot_mb03vd(n, K, 1, n, &z[0], n, n, &dwork[mem_base], n-1, &dwork[0]);
    t = z;

    slicot_mb03vy(n, K, 1, n, &z[0], n, n, &dwork[mem_base], n-1, &dwork[0], mem_needed);

    // Set numerical zeros to zero
    if (num_zero>0) {
      for (int k = 0;k<t.size();++k) {
        double &r = t[k];
        if (fabs(r)<num_zero) r = 0.0;
      }
    }

    slicot_mb03wd('S', 'V', n, K, 1, n, 1, n, &t[0], n, n, &z[0], n, n,
                  &eig_real[0], &eig_imag[0], &dwork[0], mem_needed);

  }

  CASADI_SLICOT_INTERFACE_EXPORT
  void slicot_periodic_schur(const std::vector< Matrix<double> > & a,
                             std::vector< Matrix<double> > & t, std::vector< Matrix<double> > & z,
                             std::vector< double > & eig_real, std::vector< double > & eig_imag,
                             double num_zero) {
    int K = a.size();
    int n = a[0].size1();
    for (int k=0;k<K;++k) {
      casadi_assert_message(a[k].isSquare(), "a must be square");
      casadi_assert_message(a[k].size1()==n, "a must be n-by-n");
      casadi_assert_message(a[k].isDense(), "a must be dense");
    }


    std::vector<double> a_data(n*n*K);
    // Copy data into consecutive structure
    for (int k=0;k<K;++k) {
      std::copy(a[k].begin(), a[k].end(), a_data.begin()+k*n*n);
    }

    std::vector<double> t_data(n*n*K);
    std::vector<double> z_data(n*n*K);

    slicot_periodic_schur(n, K, a_data, t_data, z_data, eig_real, eig_imag, num_zero);

    t.resize(K);
    z.resize(K);
    for (int k=0;k<K;++k) {
      t[k] = DMatrix::zeros(n, n);
      std::copy(t_data.begin()+k*n*n, t_data.begin()+(k+1)*n*n, t[k].begin());
      z[k] = DMatrix::zeros(n, n);
      std::copy(z_data.begin()+k*n*n, z_data.begin()+(k+1)*n*n, z[k].begin());
    }
  }

} // namespace casadi

/// \endcond
