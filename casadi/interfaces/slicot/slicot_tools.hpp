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


#ifndef CASADI_SLICOT_TOOLS_HPP
#define CASADI_SLICOT_TOOLS_HPP

#include "../../core/matrix/matrix.hpp"
#include <casadi/interfaces/slicot/casadi_slicot_interface_export.h>

namespace casadi {

/// \cond INTERNAL
#ifndef SWIG
void slicot_mb03vd(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2, double * tau,
                   int ldtau, double * dwork=0);

void slicot_mb03vy(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2,
                   const double * tau, int ldtau, double * dwork=0, int ldwork=0);

void slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi, int iloz, int ihiz,
                   double *h, int ldh1, int ldh2, double* z, int ldz1, int ldz2, double* wr,
                   double *wi, double * dwork=0, int ldwork=0);


void slicot_periodic_schur(int n, int K, const std::vector< double > & a,
                           std::vector< double > & t, std::vector< double > & z,
                           std::vector<double> &eig_real, std::vector<double> &eig_imag,
                           double num_zero=0);

CASADI_SLICOT_INTERFACE_EXPORT void slicot_periodic_schur(
                           int n, int K, const std::vector< double > & a, std::vector< double > & t,
                           std::vector< double > & z, std::vector<double> &dwork,
                           std::vector<double> &eig_real, std::vector<double> &eig_imag,
                           double num_zero=0);
#endif // SWIG
/// \endcond

/** \brief Obtain Periodic Schur Form of a set of matrices
*
*  Finds Z_i such that

\verbatim
          Z_1' * A_1 * Z_2 = T_1,
          Z_2' * A_2 * Z_3 = T_2,
                 ...
          Z_K' * A_K * Z_1 = T_K,
\endverbatim
*
*  with T_1 in Hessenberg form (upper triangular + one band below the diagonal)
*   and T_2..T_K  upper diagonal
*
*  with <tt>Z_k Z_k' = eye(n) = Z_k' Z_k</tt>
*
*/
CASADI_SLICOT_INTERFACE_EXPORT void slicot_periodic_schur(
    const std::vector< Matrix<double> > & A,
    std::vector< Matrix<double> > & SWIG_OUTPUT(T),
    std::vector< Matrix<double> > & SWIG_OUTPUT(Z),
    std::vector<double> &SWIG_OUTPUT(eig_real),
    std::vector<double> &SWIG_OUTPUT(eig_imag),
    double num_zero=0);

} // namespace casadi

#endif // CASADI_SLICOT_TOOLS_HPP
