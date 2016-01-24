/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file src/LAPACKReplacement.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2015
 *
 *  LAPACK replacement routines.
 */


#include <qpOASES/Utils.hpp>


extern "C" void dpotrf_(	const char *uplo, const unsigned long *_n, double *a,
							const unsigned long *_lda, long *info
							)
{
	double sum;
	long i, j, k;
	long n = (long)(*_n);
	long lda = (long)(*_lda);

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = REFER_NAMESPACE_QPOASES getSqrt( sum );
		else
		{
			a[0] = sum; /* tunnel negative diagonal element to caller */
			if (info != 0)
				*info = (long)i+1;
			return;
		}

		for( j=(i+1); j<n; ++j )
		{
			sum = a[j*lda + i];

			for( k=(i-1); k>=0; --k )
				sum -= a[k+lda*i] * a[k+lda*j];

			a[i+lda*j] = sum / a[i+lda*i];
		}
	}
	if (info != 0)
		*info = 0;
}


extern "C" void spotrf_(	const char *uplo, const unsigned long *_n, float *a,
							const unsigned long *_lda, long *info
							)
{
	float sum;
	long i, j, k;
	long n = (long)(*_n);
	long lda = (long)(*_lda);

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = (float)(REFER_NAMESPACE_QPOASES getSqrt( sum ));
		else
		{
			a[0] = sum; /* tunnel negative diagonal element to caller */
			if (info != 0)
				*info = (long)i+1;
			return;
		}

		for( j=(i+1); j<n; ++j )
		{
			sum = a[j*lda + i];

			for( k=(i-1); k>=0; --k )
				sum -= a[k+lda*i] * a[k+lda*j];

			a[i+lda*j] = sum / a[i+lda*i];
		}
	}
	if (info != 0)
		*info = 0;
}

extern "C" void dtrtrs_(	const char *UPLO, const char *TRANS, const char *DIAG,
							const unsigned long *N, const unsigned long *NRHS,
							double *A, const unsigned long *LDA, double *B, const unsigned long *LDB, long *INFO
							)
{
	; /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}

extern "C" void strtrs_(	const char *UPLO, const char *TRANS, const char *DIAG,
							const unsigned long *N, const unsigned long *NRHS,
							float *A, const unsigned long *LDA, float *B, const unsigned long *LDB, long *INFO
							)
{
	; /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}

extern "C" void dtrcon_(	const char *NORM, const char *UPLO, const char *DIAG,
							const unsigned long *N, double *A, const unsigned long *LDA,
							double *RCOND, double *WORK, const unsigned long *IWORK, long *INFO
							)
{
	; /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}

extern "C" void strcon_(	const char *NORM, const char *UPLO, const char *DIAG,
							const unsigned long *N, float *A, const unsigned long *LDA,
							float *RCOND, float *WORK, const unsigned long *IWORK, long *INFO
							)
{
	; /* Dummy. If SQProblemSchur is to be used, system LAPACK must be used */
}
