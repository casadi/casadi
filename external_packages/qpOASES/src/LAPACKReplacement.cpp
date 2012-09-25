/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2012 by Hans Joachim Ferreau, Andreas Potschka,
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
 *	\version 3.0beta
 *	\date 2007-2012
 *
 *  LAPACK replacement routines.
 */


#include <math.h>


extern "C" void dpotrf_ (const char *uplo, const unsigned long *_n, double *a,
						 const unsigned long *_lda, long *info)
{
	unsigned long n = *_n, lda = *_lda;
	double sum;
	unsigned long i, j;
	int k;

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = sqrt( sum );
		else
		{
			if (info != 0)
				*info = i+1;
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


extern "C" void spotrf_ (const char *uplo, const unsigned long *_n, float *a,
						 const unsigned long *_lda, long *info)
{
	unsigned long n = *_n, lda = *_lda;
	float sum;
	unsigned long i, j;
	int k;

	for( i=0; i<n; ++i )
	{
		/* j == i */
		sum = a[i + lda*i];

		for( k=(i-1); k>=0; --k )
			sum -= a[k+lda*i] * a[k+lda*i];

		if ( sum > 0.0 )
			a[i+lda*i] = sqrt( sum );
		else
		{
			if (info != 0)
				*info = i+1;
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

