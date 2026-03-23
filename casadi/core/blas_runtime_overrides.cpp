//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

#ifdef WITH_BLAS

#include <cstring>

#ifdef WITH_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#include "function.hpp"

namespace casadi {

  void casadi_axpy(casadi_int n, float alpha, const float* x, float* y) {
    if (!x || !y) return;
    cblas_saxpy(n, alpha, x, 1, y, 1);
  }

  void casadi_axpy(casadi_int n, double alpha, const double* x, double* y) {
    if (!x || !y) return;
    cblas_daxpy(n, alpha, x, 1, y, 1);
  }

  void casadi_copy(const float* x, casadi_int n, float* y) {
    if (y) {
      if (x) {
        cblas_scopy(n, x, 1, y, 1);
      } else {
        memset(y, 0, n * sizeof(float));
      }
    }
  }

  void casadi_copy(const double* x, casadi_int n, double* y) {
    if (y) {
      if (x) {
        cblas_dcopy(n, x, 1, y, 1);
      } else {
        memset(y, 0, n * sizeof(double));
      }
    }
  }

  float casadi_dot(casadi_int n, const float* x, const float* y) {
    return cblas_sdot(n, x, 1, y, 1);
  }

  double casadi_dot(casadi_int n, const double* x, const double* y) {
    return cblas_ddot(n, x, 1, y, 1);
  }

  float casadi_norm_1(casadi_int n, const float* x) {
    return cblas_sasum(n, x, 1);
  }

  double casadi_norm_1(casadi_int n, const double* x) {
    return cblas_dasum(n, x, 1);
  }

  float casadi_norm_2(casadi_int n, const float* x) {
    return cblas_snrm2(n, x, 1);
  }

  double casadi_norm_2(casadi_int n, const double* x) {
    return cblas_dnrm2(n, x, 1);
  }

  float casadi_norm_inf(casadi_int n, const float* x) {
    return cblas_isamax(n, x, 1);
  }

  double casadi_norm_inf(casadi_int n, const double* x) {
    return cblas_idamax(n, x, 1);
  }

  void casadi_scal(casadi_int n, float alpha, float* x) {
    if (!x) return;
    cblas_sscal(n, alpha, x, 1);
  }

  void casadi_scal(casadi_int n, double alpha, double* x) {
    if (!x) return;
    cblas_dscal(n, alpha, x, 1);
  }

}

#endif