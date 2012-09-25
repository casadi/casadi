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
 *	\file src/Utils.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches, Eckhard Arnold
 *	\version 3.0beta
 *	\date 2007-2012
 *
 *	Implementation of some utilities for working with the different QProblem classes.
 */


#include <math.h>

#if defined(__WIN32__) || defined(WIN32)
  #include <windows.h>
#elif defined(LINUX)
  #include <sys/stat.h>
  #include <sys/time.h>
#endif

#ifdef __MATLAB__
  #include "mex.h"
#endif


#include <qpOASES/Utils.hpp>


BEGIN_NAMESPACE_QPOASES


/*
 *	p r i n t
 */
returnValue print( const real_t* const v, int n )
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	int i;
	char myPrintfString[160];

	/* Print a vector. */
	myPrintf( "[\t" );
	for( i=0; i<n; ++i )
	{
		snprintf( myPrintfString,160," %.16e\t", v[i] );
		myPrintf( myPrintfString );
	}
	myPrintf( "]\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue print(	const real_t* const v, int n,
					const int* const V_idx
					)
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	int i;
	char myPrintfString[160];

	/* Print a permuted vector. */
	myPrintf( "[\t" );
	for( i=0; i<n; ++i )
	{
		snprintf( myPrintfString,160," %.16e\t", v[ V_idx[i] ] );
		myPrintf( myPrintfString );
	}
	myPrintf( "]\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue print(	const real_t* const v, int n,
					const char* name
					)
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	char myPrintfString[160];

	/* Print vector name ... */
	snprintf( myPrintfString,160,"%s = ", name );
	myPrintf( myPrintfString );
	#endif
	#endif

	/* ... and the vector itself. */
	return print( v, n );
}

/*
 *	p r i n t
 */
returnValue print( const real_t* const M, int nrow, int ncol )
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	int i;

	/* Print a matrix as a collection of row vectors. */
	for( i=0; i<nrow; ++i )
		print( &(M[i*ncol]), ncol );
	myPrintf( "\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue print(	const real_t* const M, int nrow, int ncol,
					const int* const ROW_idx, const int* const COL_idx
					)
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	int i;

	/* Print a permuted matrix as a collection of permuted row vectors. */
	for( i=0; i<nrow; ++i )
		print( &( M[ ROW_idx[i]*ncol ] ), ncol, COL_idx );
	myPrintf( "\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue print(	const real_t* const M, int nrow, int ncol,
					const char* name
					)
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	char myPrintfString[160];

	/* Print matrix name ... */
	snprintf( myPrintfString,160,"%s = ", name );
	myPrintf( myPrintfString );
	#endif
	#endif

	/* ... and the matrix itself. */
	return print( M, nrow, ncol );
}


/*
 *	p r i n t
 */
returnValue print( const int* const index, int n )
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	int i;
	char myPrintfString[160];

	/* Print a indexlist. */
	myPrintf( "[\t" );
	for( i=0; i<n; ++i )
	{
		snprintf( myPrintfString,160," %d\t", index[i] );
		myPrintf( myPrintfString );
	}
	myPrintf( "]\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue print(	const int* const index, int n,
					const char* name
					)
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	char myPrintfString[160];

	/* Print indexlist name ... */
	snprintf( myPrintfString,160,"%s = ", name );
	myPrintf( myPrintfString );
	#endif
	#endif

	/* ... and the indexlist itself. */
	return print( index, n );
}


/*
 *	m y P r i n t f
 */
returnValue myPrintf( const char* s )
{
	#ifndef __SUPPRESSANYOUTPUT__
	#ifndef __XPCTARGET__
	  #ifdef __MATLAB__
	  mexPrintf( s );
	  #else
	  FILE* outputfile = getGlobalMessageHandler( )->getOutputFile( );
	  if ( outputfile == 0 )
		  return THROWERROR( RET_NO_GLOBAL_MESSAGE_OUTPUTFILE );

	  fprintf( outputfile, "%s", s );
	  #endif
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t C o p y r i g h t N o t i c e
 */
returnValue printCopyrightNotice( )
{
	#ifndef __SUPPRESSANYOUTPUT__
		#ifndef __XPCTARGET__
		#ifndef __DSPACE__
		#ifndef __NO_COPYRIGHT__
		myPrintf( "\nqpOASES -- An Implementation of the Online Active Set Strategy.\nCopyright (C) 2007-2012 by Hans Joachim Ferreau, Andreas Potschka,\nChristian Kirches et al. All rights reserved.\n\nqpOASES is distributed under the terms of the \nGNU Lesser General Public License 2.1 in the hope that it will be \nuseful, but WITHOUT ANY WARRANTY; without even the implied warranty \nof MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. \nSee the GNU Lesser General Public License for more details.\n\n" );
		#endif
		#endif
		#endif
	#endif
	return SUCCESSFUL_RETURN;
}


/*
 *	r e a d F r o m F i l e
 */
returnValue readFromFile(	real_t* data, int nrow, int ncol,
							const char* datafilename
							)
{
	#ifndef __XPCTARGET__
	int i, j;
	real_t float_data;
	FILE* datafile;

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		char errstr[80];
		snprintf( errstr,80,"(%s)",datafilename );
		return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_OPEN_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
	}

	/* 2) Read data from file. */
	for( i=0; i<nrow; ++i )
	{
		for( j=0; j<ncol; ++j )
		{
			#ifdef __USE_SINGLE_PRECISION__
			if ( fscanf( datafile, "%f ", &float_data ) == 0 )
			#else
			if ( fscanf( datafile, "%lf ", &float_data ) == 0 )
			#endif /* __USE_SINGLE_PRECISION__ */
			{
				fclose( datafile );
				char errstr[80];
				snprintf( errstr,80,"(%s)",datafilename );
				return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_READ_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			}
			data[i*ncol + j] = ( (real_t) float_data );
		}
	}

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;
	#else

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	r e a d F r o m F i l e
 */
returnValue readFromFile(	real_t* data, int n,
							const char* datafilename
							)
{
	return readFromFile( data, n, 1, datafilename );
}



/*
 *	r e a d F r o m F i l e
 */
returnValue readFromFile(	int* data, int n,
							const char* datafilename
							)
{
	#ifndef __XPCTARGET__
	int i;
	FILE* datafile;

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		char errstr[80];
		snprintf( errstr,80,"(%s)",datafilename );
		return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_OPEN_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
	}

	/* 2) Read data from file. */
	for( i=0; i<n; ++i )
	{
		if ( fscanf( datafile, "%d\n", &(data[i]) ) == 0 )
		{
			fclose( datafile );
			char errstr[80];
			snprintf( errstr,80,"(%s)",datafilename );
			return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_READ_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;
	#else

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	w r i t e I n t o F i l e
 */
returnValue writeIntoFile(	const real_t* const data, int nrow, int ncol,
							const char* datafilename, BooleanType append
							)
{
	#ifndef __XPCTARGET__
	int i, j;
	FILE* datafile;

	/* 1) Open file. */
	if ( append == BT_TRUE )
	{
		/* append data */
		if ( ( datafile = fopen( datafilename, "a" ) ) == 0 )
		{
			char errstr[80];
			snprintf( errstr,80,"(%s)",datafilename );
			return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_OPEN_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}
	else
	{
		/* do not append data */
		if ( ( datafile = fopen( datafilename, "w" ) ) == 0 )
		{
			char errstr[80];
			snprintf( errstr,80,"(%s)",datafilename );
			return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_OPEN_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}

	/* 2) Write data into file. */
	for( i=0; i<nrow; ++i )
	{
		for( j=0; j<ncol; ++j )
		 	fprintf( datafile, "%.16e ", data[i*ncol+j] );

		fprintf( datafile, "\n" );
	}

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;
	#else

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	w r i t e I n t o F i l e
 */
returnValue writeIntoFile(	const real_t* const data, int n,
							const char* datafilename, BooleanType append
							)
{
	return writeIntoFile( data,1,n,datafilename,append );
}


/*
 *	w r i t e I n t o F i l e
 */
returnValue writeIntoFile(	const int* const integer, int n,
							const char* datafilename, BooleanType append
							)
{
	#ifndef __XPCTARGET__
	int i;

	FILE* datafile;

	/* 1) Open file. */
	if ( append == BT_TRUE )
	{
		/* append data */
		if ( ( datafile = fopen( datafilename, "a" ) ) == 0 )
		{
			char errstr[80];
			snprintf( errstr,80,"(%s)",datafilename );
			return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_OPEN_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}
	else
	{
		/* do not append data */
		if ( ( datafile = fopen( datafilename, "w" ) ) == 0 )
		{
			char errstr[80];
			snprintf( errstr,80,"(%s)",datafilename );
			return getGlobalMessageHandler( )->throwError( RET_UNABLE_TO_OPEN_FILE,errstr,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}

	/* 2) Write data into file. */
	for( i=0; i<n; ++i )
		fprintf( datafile, "%d\n", integer[i] );

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;
	#else

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	g e t C P U t i m e
 */
real_t getCPUtime( )
{
	real_t current_time = -1.0;

	#if defined(__WIN32__) || defined(WIN32)
	LARGE_INTEGER counter, frequency;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&counter);
	current_time = ((real_t) counter.QuadPart) / ((real_t) frequency.QuadPart);
	#elif defined(LINUX)
	struct timeval theclock;
	gettimeofday( &theclock,0 );
	current_time = 1.0*theclock.tv_sec + 1.0e-6*theclock.tv_usec;
	#endif

	return current_time;
}


/*
 *	g e t N o r m
 */
real_t getNorm( const real_t* const v, int n )
{
	int i;

	real_t norm = 0.0;

	for( i=0; i<n; ++i )
		norm += v[i]*v[i];

	return sqrt( norm );
}



/*
 *	g e t K K T R e s i d u a l
 */
void getKKTResidual(	int nV, int nC,
						const real_t* const H, const real_t* const g,
						const real_t* const A, const real_t* const lb, const real_t* const ub,
						const real_t* const lbA, const real_t* const ubA,
						const real_t* const x, const real_t* const y,
						real_t& stat, real_t& feas, real_t& cmpl)
{
	/* Tolerance for dual variables considered zero. */
	const real_t dualActiveTolerance = 1.0e3 * EPS;

	int i, j;
	real_t sum, prod;

	/* Initialize residuals */
	stat = feas = cmpl = 0.0;

	/* Pointers to data and solution of current QP */
	const real_t* gCur   = g;
	const real_t* lbCur  = lb;
	const real_t* ubCur  = ub;
	const real_t* lbACur = lbA ? lbA : 0;
	const real_t* ubACur = ubA ? ubA : 0;
	const real_t* xCur   = x;
	const real_t* yCur   = y;

	/* check stationarity */
	for (i = 0; i < nV; i++)
	{
		/* g term and variable bounds dual term */
		sum = gCur[i] - yCur[i];

		/* H*x term */
		for (j = 0; j < nV; j++) sum += H[i*nV+j] * xCur[j];
		/* A'*y term */
		for (j = 0; j < nC; j++) sum -= A[j*nV+i] * yCur[j+nV];

		/* update stat */
		if (fabs(sum) > stat) stat = fabs(sum);
	}

	/* check primal feasibility and complementarity */
	/* variable bounds */
	for (i = 0; i < nV; i++)
	{
		/* feasibility */
		if (lbCur[i] - xCur[i] > feas) 
			feas = lbCur[i] - xCur[i];
		if (xCur[i] - ubCur[i] > feas) 
			feas = xCur[i] - ubCur[i];

		/* complementarity */
		prod = 0.0;
		if (yCur[i] > dualActiveTolerance) /* lower bound */
			prod = (xCur[i] - lbCur[i]) * yCur[i];
		if (yCur[i] < -dualActiveTolerance) /* upper bound */
			prod = (xCur[i] - ubCur[i]) * yCur[i];
		if (fabs(prod) > cmpl) cmpl = fabs(prod);
	}
	/* A*x bounds */
	for (i = 0; i < nC; i++)
	{
		/* compute sum = (A*x)_i */
		sum = 0.0;
		for (j = 0; j < nV; j++) 
			sum += A[i*nV+j] * xCur[j];

		/* feasibility */
		if (lbACur[i] - sum > feas) 
			feas = lbACur[i] - sum;
		if (sum - ubACur[i] > feas) 
			feas = sum - ubACur[i];

		/* complementarity */
		prod = 0.0;
		if (yCur[nV+i] > dualActiveTolerance) /* lower bound */
			prod = (sum - lbACur[i]) * yCur[nV+i];
		if (yCur[nV+i] < -dualActiveTolerance) /* upper bound */
			prod = (sum - ubACur[i]) * yCur[nV+i];
		if (fabs(prod) > cmpl) cmpl = fabs(prod);
	}
}


returnValue convertBooleanTypeToString( BooleanType value, char* const string )
{
	if ( value == BT_TRUE )
		snprintf( string,20,"BT_TRUE" );
	else
		snprintf( string,20,"BT_FALSE" );

	return SUCCESSFUL_RETURN;
}


returnValue convertSubjectToStatusToString( SubjectToStatus value, char* const string )
{
	switch( value )
	{
		case ST_INACTIVE:
			snprintf( string,20,"ST_INACTIVE" );
			break;

		case ST_LOWER:
			snprintf( string,20,"ST_LOWER" );
			break;

		case ST_UPPER:
			snprintf( string,20,"ST_UPPER" );
			break;

		case ST_UNDEFINED:
			snprintf( string,20,"ST_UNDEFINED" );
			break;
			
		case ST_INFEASIBLE_LOWER:
			snprintf( string,20,"ST_INFEASIBLE_LOWER" );
			break;

		case ST_INFEASIBLE_UPPER:
			snprintf( string,20,"ST_INFEASIBLE_UPPER" );
			break;
	}

	return SUCCESSFUL_RETURN;
}


returnValue convertPrintLevelToString( PrintLevel value, char* const string )
{
	switch( value )
	{
		case PL_NONE:
			snprintf( string,20,"PL_NONE" );
			break;

		case PL_LOW:
			snprintf( string,20,"PL_LOW" );
			break;

		case PL_MEDIUM:
			snprintf( string,20,"PL_MEDIUM" );
			break;

		case PL_HIGH:
			snprintf( string,20,"PL_HIGH" );
			break;
			
		case PL_TABULAR:
			snprintf( string,20,"PL_TABULAR" );
			break;
	}

	return SUCCESSFUL_RETURN;
}



#ifdef __DEBUG__
extern "C" void gdb_printmat(const char *fname, real_t *M, int n, int m, int ldim)
{
	int i, j;
	FILE *fid;

	fid = fopen(fname, "wt");
	if (!fid) 
	{
		perror("Error opening file: ");
		return;
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
			fprintf(fid, " %23.16e", M[j*ldim+i]);
		fprintf(fid, "\n");
	}
	fclose(fid);
}
#endif /* __DEBUG__ */



END_NAMESPACE_QPOASES


/*
 *	end of file
 */
