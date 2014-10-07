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


#include "matrix_tools.hpp"
#include "../std_vector_tools.hpp"
#include "../function/linear_solver.hpp"

using namespace std;

namespace casadi {

  Matrix<double> solve(const Matrix<double>& A, const Matrix<double>& b,
                       const std::string& lsolver, const Dictionary& dict) {
    LinearSolver mysolver(lsolver, A.sparsity(), b.size2());
    mysolver.setOption(dict);
    mysolver.init();
    mysolver.setInput(A, LINSOL_A);
    mysolver.setInput(b, LINSOL_B);
    mysolver.prepare();
    mysolver.solve(false);

    return mysolver.output(LINSOL_X);
  }


  Matrix<double> pinv(const Matrix<double>& A, const std::string& lsolver,
                      const Dictionary& dict) {
    if (A.size1()>=A.size2()) {
      return solve(mul(A.T(), A), A.T(), lsolver, dict);
    } else {
      return solve(mul(A, A.T()), A, lsolver, dict).T();
    }
  }

  double norm_inf_mul_nn(const Matrix<double> &B,
                         const Matrix<double> &A,
                         std::vector<double>& Dwork,
                         std::vector<int>& Iwork) {

    // Note: because the algorithm works with compressed row storage,
    // we have x=B and y=A
    double res = 0;

    casadi_assert_message(A.size1()==B.size2(), "Dimension error. Got " << B.dimString()
                          << " times " << A.dimString() << ".");


    int n_row = A.size2();
    int n_col = B.size1();

    casadi_assert_message(Dwork.size()>=n_col,
      "We need a bigger work vector (>=" << n_col << "), but got "<< Dwork.size() <<".");
    casadi_assert_message(Iwork.size()>=n_row+1+n_col,
      "We need a bigger work vector (>=" << n_row+1+n_col << "), but got "<< Iwork.size() <<".");

    const std::vector<int> &Aj = A.row();
    const std::vector<int> &Ap = A.colind();
    const std::vector<double> &Ax = A.data();

    const std::vector<int> &Bj = B.row();
    const std::vector<int> &Bp = B.colind();
    const std::vector<double> &Bx = B.data();

    int *Cp = &Iwork[0];
    int *mask = &Iwork[n_row+1];

    // Implementation borrowed from Scipy's sparsetools/csr.h

    // Pass 1

    // method that uses O(n) temp storage
    std::fill(mask, mask+n_col, -1);

    Cp[0] = 0;
    int nnz = 0;

    for (int i = 0; i < n_row; i++) {
      int row_nnz = 0;
      for (int jj = Ap[i]; jj < Ap[i+1]; jj++) {
        int j = Aj[jj];
        for (int kk = Bp[j]; kk < Bp[j+1]; kk++) {
          int k = Bj[kk];
          if (mask[k] != i) {
            mask[k] = i;
            row_nnz++;
          }
        }
      }
      int next_nnz = nnz + row_nnz;

      nnz = next_nnz;
      Cp[i+1] = nnz;
    }

    // Pass 2
    int *next = &Iwork[n_row+1];
    std::fill(next, next+n_col, -1);

    double* sums = &Dwork[0];
    std::fill(sums, sums+n_col, 0);

    nnz = 0;

    Cp[0] = 0;

    for (int i = 0; i < n_row; i++) {
        int head   = -2;
        int length =  0;

        int jj_start = Ap[i];
        int jj_end   = Ap[i+1];
        for (int jj = jj_start; jj < jj_end; jj++) {
            int j = Aj[jj];
            double v = Ax[jj];

            int kk_start = Bp[j];
            int kk_end   = Bp[j+1];
            for (int kk = kk_start; kk < kk_end; kk++) {
                int k = Bj[kk];

                sums[k] += v*Bx[kk];

                if (next[k] == -1) {
                    next[k] = head;
                    head  = k;
                    length++;
                }
            }
        }

        for (int jj = 0; jj < length; jj++) {

            if (sums[head] != 0) {
                res = std::max(res, abs(sums[head]));
                nnz++;
            }

            int temp = head;
            head = next[head];

            next[temp] = -1; //clear arrays
            sums[temp] =  0;
        }

        Cp[i+1] = nnz;
    }

    return res;

  }


} // namespace casadi


