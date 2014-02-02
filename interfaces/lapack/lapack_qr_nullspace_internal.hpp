/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef LAPACK_QR_NULLSPACE_INTERNAL_HPP
#define LAPACK_QR_NULLSPACE_INTERNAL_HPP

#include "lapack_qr_nullspace.hpp"
#include "../../symbolic/fx/nullspace_internal.hpp"

namespace CasADi{

  /** @copydoc Nullspace_doc
      \author Joris Gillis 
      \date 2014
  */
  class LapackQRNullspaceInternal : public NullspaceInternal{
  public:
  
    /** \brief  Constructor */
    LapackQRNullspaceInternal(const CRSSparsity & A_sp, int nfwd=0, int nadj=0);
  
    /** \brief  Destructor */
    virtual ~LapackQRNullspaceInternal();
  
    /** \brief  Clone */
    virtual LapackQRNullspaceInternal* clone() const;
    
    /** \brief  Create a new solver */
    virtual LapackQRNullspaceInternal* create(const CRSSparsity & A_sp_) const{ return new LapackQRNullspaceInternal(A_sp_);}

    /** \brief  initialize */
    virtual void init();

    /** \brief  Integrate */
    virtual void evaluate();
    
    /// Generate a function that calculates nfwd forward derivatives and nadj adjoint derivatives
    virtual FX getDerivative(int nfwd, int nadj);
    
  private:
  
    /// Number of fwd directions
    int nfwd_;
    
    /// Number of fwd directions
    int nadj_;
    
    std::vector<double> tempA_;
    std::vector<double> tempB_;
    std::vector<double> tempC_;
  
    std::vector<double> tau_;
    
    std::vector<double> work_;
    std::vector<double> work2_;
    
    Matrix<double> sselector_;
    
  };
  
    /// QR-Factorize dense matrix (lapack)
  extern "C" void dgeqrf_(const int *m,const int *n, double *a,const int *lda, double *tau, double *work,const int *lwork, int *info);
  
  extern "C" void dormqr_(const char *side,const  char*trans,const  int *m,const  int *n,const  int *K,const  double *a,const  int *lda,const  double *tau, double* C,const  int *ldc, double *work,const  int *lwork, int *info);
  
  extern "C" void dtrtrs_(const char *uplo, const char* trans, const char* diag, const int *n, const int *nrhs, double *a, const int *lda, const double *b, int *ldb, int *info);

  
} // namespace CasADi

#endif // LAPACK_QR_NULLSPACE_INTERNAL_HPP
