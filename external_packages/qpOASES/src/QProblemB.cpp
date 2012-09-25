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
 *	\file src/QProblemB.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0beta
 *	\date 2007-2012
 *
 *	Implementation of the QProblemB class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 */


#include <qpOASES/QProblemB.hpp>

BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	Q P r o b l e m B
 */
QProblemB::QProblemB( )
{
	/* print copyright notice */
	if (options.printLevel > PL_NONE)
		printCopyrightNotice( );

	/* reset global message handler */
	getGlobalMessageHandler( )->reset( );

	freeHessian = BT_FALSE;
	H = 0;

	g = 0;
	lb = 0;
	ub = 0;

	R = 0;
	haveCholesky = BT_FALSE;

	x = 0;
	y = 0;

	tau = 0.0;

	hessianType = HST_UNKNOWN;
	isRegularised = BT_FALSE;

	infeasible  = BT_FALSE;
	unbounded   = BT_FALSE;

	status = QPS_NOTINITIALISED;

	count = 0;

	ramp0 = options.initialRamping;
	ramp1 = options.finalRamping;
	rampOffset = 0;

	delta_xFR_TMP = 0;

	setPrintLevel( options.printLevel );
}


/*
 *	Q P r o b l e m B
 */
QProblemB::QProblemB( int _nV, HessianType _hessianType )
{
	/* print copyright notice */
	if (options.printLevel > PL_NONE)
		printCopyrightNotice( );

	/* consistency check */
	if ( _nV <= 0 )
	{
		_nV = 1;
		THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* reset global message handler */
	getGlobalMessageHandler( )->reset( );

	freeHessian = BT_FALSE;
	H = 0;

	g = new real_t[_nV];
	lb = new real_t[_nV];
	ub = new real_t[_nV];

	bounds.init( _nV );

	R = new real_t[_nV*_nV];

	x = new real_t[_nV];
	y = new real_t[_nV];

	tau = 0.0;

	hessianType = _hessianType;
	isRegularised = BT_FALSE;

	infeasible  = BT_FALSE;
	unbounded   = BT_FALSE;

	status = QPS_NOTINITIALISED;

	count = 0;

	ramp0 = options.initialRamping;
	ramp1 = options.finalRamping;
	rampOffset = 0;

	delta_xFR_TMP = new real_t[_nV];

	setPrintLevel( options.printLevel );

	flipper.init( _nV );
}


/*
 *	Q P r o b l e m B
 */
QProblemB::QProblemB( const QProblemB& rhs )
{
	freeHessian = BT_FALSE;
	H = 0;

	copy( rhs );
}


/*
 *	~ Q P r o b l e m B
 */
QProblemB::~QProblemB( )
{
	clear( );

	/* reset global message handler */
	getGlobalMessageHandler( )->reset( );
}


/*
 *	o p e r a t o r =
 */
QProblemB& QProblemB::operator=( const QProblemB& rhs )
{
	if ( this != &rhs )
	{
		clear( );
		copy( rhs );
	}

	return *this;
}


/*
 *	r e s e t
 */
returnValue QProblemB::reset( )
{
	int i;
	int nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Reset bounds. */
	bounds.init( nV );

	/* 2) Reset Cholesky decomposition. */
	for( i=0; i<nV*nV; ++i )
		R[i] = 0.0;
	
	haveCholesky = BT_FALSE;

	/* 3) Reset steplength and status flags. */
	tau = 0.0;

	hessianType = HST_UNKNOWN;
	isRegularised = BT_FALSE;

	infeasible  = BT_FALSE;
	unbounded   = BT_FALSE;

	status = QPS_NOTINITIALISED;

	ramp0 = options.initialRamping;
	ramp1 = options.finalRamping;
	rampOffset = 0;

	return SUCCESSFUL_RETURN;
}


/*
 *	i n i t
 */
returnValue QProblemB::init(	SymmetricMatrix *_H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int& nWSR, real_t* const cputime
								)
{
	if ( getNV( ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine (without any additional information). */
	return solveInitialQP( 0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB::init(	const real_t* const _H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int& nWSR, real_t* const cputime
								)
{
	if ( getNV( ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine (without any additional information). */
	return solveInitialQP( 0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB::init(	const char* const H_file, const char* const g_file,
								const char* const lb_file, const char* const ub_file,
								int& nWSR, real_t* const cputime
								)
{
	if ( getNV( ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency check. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	/* 2) Setup QP data from files. */
	if ( setupQPdataFromFile( H_file,g_file,lb_file,ub_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* 3) Call to main initialisation routine (without any additional information). */
	return solveInitialQP( 0,0,0, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB::init( 	SymmetricMatrix *_H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds
								)
{
	int i;
	int nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( guessedBounds->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	/* exclude this possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return solveInitialQP( xOpt,yOpt,guessedBounds, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB::init( 	const real_t* const _H, const real_t* const _g,
								const real_t* const _lb, const real_t* const _ub,
								int& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds
								)
{
	int i;
	int nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	if ( guessedBounds != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( guessedBounds->getStatus( i ) == ST_UNDEFINED )
				return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	/* exclude this possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data. */
	if ( setupQPdata( _H,_g,_lb,_ub ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 3) Call to main initialisation routine. */
	return solveInitialQP( xOpt,yOpt,guessedBounds, nWSR,cputime );
}


/*
 *	i n i t
 */
returnValue QProblemB::init( 	const char* const H_file, const char* const g_file,
								const char* const lb_file, const char* const ub_file,
								int& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds
								)
{
	int i;
	int nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* 1) Consistency checks. */
	if ( isInitialised( ) == BT_TRUE )
	{
		THROWWARNING( RET_QP_ALREADY_INITIALISED );
		reset( );
	}

	for( i=0; i<nV; ++i )
	{
		if ( guessedBounds->getStatus( i ) == ST_UNDEFINED )
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* exclude this possibility in order to avoid inconsistencies */
	if ( ( xOpt == 0 ) && ( yOpt != 0 ) && ( guessedBounds != 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* 2) Setup QP data from files. */
	if ( setupQPdataFromFile( H_file,g_file,lb_file,ub_file ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* 3) Call to main initialisation routine. */
	return solveInitialQP( xOpt,yOpt,guessedBounds, nWSR,cputime );
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB::hotstart(	const real_t* const g_new,
									const real_t* const lb_new, const real_t* const ub_new,
									int& nWSR, real_t* const cputime
									)
{
	if ( getNV( ) == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	++count;

	returnValue returnvalue = SUCCESSFUL_RETURN;
	int nFar, i;
	int nV = getNV ();

	int nWSR_max = nWSR;
	int nWSR_performed = 0;

	real_t cputime_remaining = INFTY;
	real_t cputime_needed = 0.0;

	real_t farbound = options.initialFarBounds;
	real_t t, rampVal;

	if ( options.enableFarBounds == BT_FALSE )
	{
		/*if (haveCholesky == BT_FALSE)
		{
			returnvalue = computeInitialCholesky();
			if (returnvalue != SUCCESSFUL_RETURN)
				return THROWERROR(returnvalue);
		}*/

		if ( usingRegularisation( ) == BT_TRUE )
			returnvalue = solveRegularisedQP( g_new,lb_new,ub_new, nWSR,cputime,0 );
		else
			returnvalue = solveQP( g_new,lb_new,ub_new, nWSR,cputime,0 );
	}
	else
	{
		real_t *ub_new_far = new real_t[nV];
		real_t *lb_new_far = new real_t[nV];

		/* possibly extend initial far bounds to largest bound/constraint data */
		if (ub_new)
			for (i = 0; i < nV; i++)
				if (ub_new[i] < INFTY && ub_new[i] > farbound) farbound = ub_new[i];
		if (lb_new)
			for (i = 0; i < nV; i++)
				if (lb_new[i] > -INFTY && lb_new[i] < -farbound) farbound = -lb_new[i];

		if ( options.enableRamping == BT_TRUE )
		{
			/* TODO: encapsule this part to avoid code duplication! */
			for ( i=0; i<nV; ++i )
			{
				t = static_cast<real_t>((i + rampOffset) % nV) / static_cast<real_t>(nV-1);
				rampVal = farbound + (1.0-t) * ramp0 + t * ramp1;

				if ( ( lb_new == 0 ) || ( lb_new[i] <= -rampVal ) )
					lb_new_far[i] = - rampVal;
				else
					lb_new_far[i] = lb_new[i];
				if ( ( ub_new == 0 ) || ( ub_new[i] >= rampVal ) )
					ub_new_far[i] = rampVal;
				else
					ub_new_far[i] = ub_new[i];
			}
		}
		else
		{
			for ( i=0; i<nV; ++i )
			{
				lb_new_far[i] = lb_new[i];
				ub_new_far[i] = ub_new[i];
			}
		}

		/*if (haveCholesky == BT_FALSE)
		{
			returnvalue = computeInitialCholesky();
			if (returnvalue != SUCCESSFUL_RETURN)
				goto farewell;
		}*/

		for ( ;; )
		{
			nWSR = nWSR_max;
			if ( cputime != 0 )
				cputime_remaining = *cputime - cputime_needed;

			if ( usingRegularisation( ) == BT_TRUE )
				returnvalue = solveRegularisedQP( g_new,lb_new_far,ub_new_far, nWSR,&cputime_remaining,nWSR_performed );
			else
				returnvalue = solveQP( g_new,lb_new_far,ub_new_far, nWSR,&cputime_remaining,nWSR_performed );

			nWSR_performed  = nWSR;
			cputime_needed += cputime_remaining;

			/* Check for active far-bounds and move them away */
			nFar = 0;
			farbound *= options.growFarBounds;

			real_t maxFarbound = 1e20;
			if ( infeasible == BT_TRUE )
			{
				if ( farbound > maxFarbound )
				{
					returnvalue = RET_HOTSTART_STOPPED_INFEASIBILITY;
					goto farewell;
				}

				if ( options.enableRamping == BT_TRUE )
				{
					/* TODO: encapsule this part to avoid code duplication! */
					for ( i=0; i<nV; ++i )
					{
						t = static_cast<real_t>((i + rampOffset) % nV) / static_cast<real_t>(nV-1);
						rampVal = farbound + (1.0-t) * ramp0 + t * ramp1;

						if ( lb_new == 0 )
							lb_new_far[i] = - rampVal;
						else
							lb_new_far[i] = getMax (- rampVal, lb_new[i]);

						if ( ub_new == 0 )
							ub_new_far[i] = rampVal;
						else
							ub_new_far[i] = getMin (rampVal, ub_new[i]);
					}
				}
				else
				{
					for ( i=0; i<nV; ++i )
					{
						lb_new_far[i] = lb_new[i];
						ub_new_far[i] = ub_new[i];
					}
				}
			}
			else if ( status == QPS_SOLVED )
			{
				real_t tol = farbound * options.boundTolerance;
				status = QPS_HOMOTOPYQPSOLVED;

				nFar = 0;
				for ( i=0; i<nV; ++i )
				{
					if ( ( ( lb_new == 0 ) || ( lb_new_far[i] > lb_new[i] ) ) && ( getAbs ( lb_new_far[i] - x[i] ) < tol ) )
						++nFar;
					if ( ( ( ub_new == 0 ) || ( ub_new_far[i] < ub_new[i] ) ) && ( getAbs ( ub_new_far[i] - x[i] ) < tol ) )
						++nFar;
				}

				if ( nFar == 0 )
					break;

				if ( farbound > maxFarbound )
				{
					unbounded = BT_TRUE;
					returnvalue = RET_HOTSTART_STOPPED_UNBOUNDEDNESS;
					goto farewell;
				}

				if ( options.enableRamping == BT_TRUE )
				{
					/* TODO: encapsule this part to avoid code duplication! */
					for ( i=0; i<nV; ++i )
					{
						t = static_cast<real_t>((i + rampOffset) % nV) / static_cast<real_t>(nV-1);
						rampVal = farbound + (1.0-t) * ramp0 + t * ramp1;

						if ( lb_new == 0 )
							lb_new_far[i] = - rampVal;
						else
							lb_new_far[i] = getMax (- rampVal, lb_new[i]);

						if ( ub_new == 0 )
							ub_new_far[i] = rampVal;
						else
							ub_new_far[i] = getMin (rampVal, ub_new[i]);
					}
				}
				else
				{
					for ( i=0; i<nV; ++i )
					{
						lb_new_far[i] = lb_new[i];
						ub_new_far[i] = ub_new[i];
					}
				}
			}
			else
			{
				/* some other error */
				break;
			}

			/* advance ramp offset to avoid Ramping cycles */
			rampOffset++;
		}

		farewell:
			if ( cputime != 0 )
				*cputime = cputime_needed;
			delete[] lb_new_far; delete[] ub_new_far;
	}

	return ( returnvalue != SUCCESSFUL_RETURN ) ? THROWERROR( returnvalue ) : returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB::hotstart(	const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									int& nWSR, real_t* const cputime
									)
{
	int nV  = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Allocate memory (if bounds exist). */
	real_t* g_new  = new real_t[nV];
	real_t* lb_new = 0;
	real_t* ub_new = 0;

	if ( lb_file != 0 )
		lb_new = new real_t[nV];
	if ( ub_file != 0 )
		ub_new = new real_t[nV];

	/* 2) Load new QP vectors from file. */
	returnValue returnvalue;
	returnvalue = loadQPvectorsFromFile(	g_file,lb_file,ub_file,
											g_new,lb_new,ub_new
											);
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		if ( ub_file != 0 )
			delete[] ub_new;
		if ( lb_file != 0 )
			delete[] lb_new;
		delete[] g_new;

		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* 3) Actually perform hotstart. */
	returnvalue = hotstart( g_new,lb_new,ub_new, nWSR,cputime );

	/* 4) Free memory. */
	if ( ub_file != 0 )
		delete[] ub_new;
	if ( lb_file != 0 )
		delete[] lb_new;
	delete[] g_new;

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB::hotstart(	const real_t* const g_new,
									const real_t* const lb_new, const real_t* const ub_new,
									int& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds
									)
{
	int nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );


	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = getCPUtime( );


	/* 1) Update working set according to guess for working set of bounds. */
	if ( guessedBounds != 0 )
	{
		if ( setupAuxiliaryQP( guessedBounds ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	else
	{
		/* create empty bounds for setting up auxiliary QP */
		Bounds emptyBounds( nV );
		if ( emptyBounds.setupAllFree( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( setupAuxiliaryQP( &emptyBounds ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}

	/* 2) Perform usual homotopy. */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= getCPUtime( ) - starttime;

	returnValue returnvalue = hotstart( g_new,lb_new,ub_new, nWSR,cputime );

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = getCPUtime( ) - starttime;

	return returnvalue;
}


/*
 *	h o t s t a r t
 */
returnValue QProblemB::hotstart(	const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									int& nWSR, real_t* const cputime,
									const Bounds* const guessedBounds
									)
{
	int nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	/* consistency check */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Allocate memory (if bounds exist). */
	real_t* g_new  = new real_t[nV];
	real_t* lb_new = 0;
	real_t* ub_new = 0;

	if ( lb_file != 0 )
		lb_new = new real_t[nV];
	if ( ub_file != 0 )
		ub_new = new real_t[nV];

	/* 2) Load new QP vectors from file. */
	returnValue returnvalue;
	returnvalue = loadQPvectorsFromFile(	g_file,lb_file,ub_file,
											g_new,lb_new,ub_new
											);
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		if ( ub_file != 0 )
			delete[] ub_new;
		if ( lb_file != 0 )
			delete[] lb_new;
		delete[] g_new;

		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	/* 3) Actually perform hotstart using initialised homotopy. */
	returnvalue = hotstart(	g_new,lb_new,ub_new, nWSR,cputime,
							guessedBounds
							);

	/* 4) Free memory. */
	if ( ub_file != 0 )
		delete[] ub_new;
	if ( lb_file != 0 )
		delete[] lb_new;
	delete[] g_new;

	return returnvalue;
}


/*
 *	g e t N Z
 */
int QProblemB::getNZ( ) const
{
	/* if no constraints are present: nZ=nFR */
	return getNFR( );
}


/*
 *	g e t O b j V a l
 */
real_t QProblemB::getObjVal( ) const
{
	real_t objVal;

	/* calculated optimal objective function value
	 * only if current QP has been solved */
	if ( ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED ) )
	{
		objVal = getObjVal( x );
	}
	else
	{
		objVal = INFTY;
	}

	return objVal;
}


/*
 *	g e t O b j V a l
 */
real_t QProblemB::getObjVal( const real_t* const _x ) const
{
	int i;
	int nV = getNV( );

	if ( nV == 0 )
		return 0.0;

	real_t objVal = 0.0;

	for( i=0; i<nV; ++i )
		objVal += _x[i]*g[i];

	switch ( hessianType )
	{
		case HST_ZERO:
			break;

		case HST_IDENTITY:
			for( i=0; i<nV; ++i )
				objVal += 0.5*_x[i]*_x[i];
			break;

		default:
			real_t *Hx = new real_t[nV];
			H->times(1, 1.0, _x, nV, 0.0, Hx, nV);
			for( i=0; i<nV; ++i )
				objVal += 0.5*_x[i]*Hx[i];
			delete[] Hx;
			break;
	}

	/* When using regularisation, the objective function value
	 * needs to be modified as follows:
	 * objVal = objVal - 0.5*_x*(Hmod-H)*_x - _x'*(gMod-g)
	 *        = objVal - 0.5*_x*eps*_x * - _x'*(-eps*_x)
	 *        = objVal + 0.5*_x*eps*_x */
	if ( usingRegularisation( ) == BT_TRUE )
	{
		for( i=0; i<nV; ++i )
			objVal += 0.5*_x[i]*options.epsRegularisation*_x[i];
	}

	return objVal;
}


/*
 *	g e t P r i m a l S o l u t i o n
 */
returnValue QProblemB::getPrimalSolution( real_t* const xOpt ) const
{
	int i;

	/* return optimal primal solution vector
	 * only if current QP has been solved */
	if ( ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED ) )
	{
		for( i=0; i<getNV( ); ++i )
			xOpt[i] = x[i];

		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}


/*
 *	g e t D u a l S o l u t i o n
 */
returnValue QProblemB::getDualSolution( real_t* const yOpt ) const
{
	int i;

	for( i=0; i<getNV( ); ++i )
		yOpt[i] = y[i];

	/* return optimal dual solution vector
	 * only if current QP has been solved */
	if ( ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED ) )
	{
		return SUCCESSFUL_RETURN;
	}
	else
	{
		return RET_QP_NOT_SOLVED;
	}
}


/*
 *	s e t P r i n t L e v e l
 */
returnValue QProblemB::setPrintLevel( PrintLevel _printLevel )
{
	#ifndef __MATLAB__
	if ( ( options.printLevel == PL_HIGH ) && ( options.printLevel != _printLevel ) )
		THROWINFO( RET_PRINTLEVEL_CHANGED );
	#endif

	options.printLevel = _printLevel;

	/* update message handler preferences */
 	switch ( options.printLevel )
 	{
		case PL_TABULAR:
 		case PL_NONE:
 			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			break;

		case PL_LOW:
			#ifndef __XPCTARGET__
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			#endif
			break;

		case PL_MEDIUM:
			#ifndef __XPCTARGET__
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			#endif
			break;

		default: /* PL_HIGH */
			#ifndef __XPCTARGET__
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_VISIBLE );
			#endif
			break;
 	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t P r o p e r t i e s
 */
returnValue QProblemB::printProperties( )
{
	#ifndef __XPCTARGET__
	#ifndef __DSPACE__
	/* Do not print properties if print level is set to none! */
	if ( options.printLevel == PL_NONE )
		return SUCCESSFUL_RETURN;

	char myPrintfString[80];

	myPrintf( "\n#################   qpOASES  --  QP PROPERTIES   #################\n" );
	myPrintf( "\n" );

	/* 1) Variables properties. */
	snprintf( myPrintfString,80,  "Number of Variables: %4.1d\n",getNV( ) );
	myPrintf( myPrintfString );

	if ( bounds.hasNoLower( ) == BT_TRUE )
			myPrintf( "Variables are not bounded from below.\n" );
		else
			myPrintf( "Variables are bounded from below.\n" );

	if ( bounds.hasNoUpper( ) == BT_TRUE )
			myPrintf( "Variables are not bounded from above.\n" );
		else
			myPrintf( "Variables are bounded from above.\n" );

	myPrintf( "\n" );


	/* 2) Further properties. */
	switch ( hessianType )
	{
		case HST_ZERO:
			myPrintf( "Hessian is zero matrix (i.e. actually an LP is solved).\n" );
			break;

		case HST_IDENTITY:
			myPrintf( "Hessian is identity matrix.\n" );
			break;

		case HST_POSDEF:
			myPrintf( "Hessian matrix is (strictly) positive definite.\n" );
			break;

		case HST_POSDEF_NULLSPACE:
			myPrintf( "Hessian matrix is positive definite on null space of active constraints.\n" );
			break;

		case HST_SEMIDEF:
			myPrintf( "Hessian matrix is positive semi-definite.\n" );
			break;

		default:
			myPrintf( "Hessian matrix has unknown type.\n" );
			break;
	}

	if ( infeasible == BT_TRUE )
		myPrintf( "QP was found to be infeasible.\n" );
	else
		myPrintf( "QP seems to be feasible.\n" );

	if ( unbounded == BT_TRUE )
		myPrintf( "QP was found to be unbounded from below.\n" );
	else
		myPrintf( "QP seems to be bounded from below.\n" );

	myPrintf( "\n" );


	/* 3) QP object properties. */
	switch ( status )
	{
		case QPS_NOTINITIALISED:
			myPrintf( "Status of QP object: freshly instantiated or reset.\n" );
			break;

		case QPS_PREPARINGAUXILIARYQP:
			myPrintf( "Status of QP object: an auxiliary QP is currently setup.\n" );
			break;

		case QPS_AUXILIARYQPSOLVED:
			myPrintf( "Status of QP object: an auxilary QP was solved.\n" );
			break;

		case QPS_PERFORMINGHOMOTOPY:
			myPrintf( "Status of QP object: a homotopy step is performed.\n" );
			break;

		case QPS_HOMOTOPYQPSOLVED:
			myPrintf( "Status of QP object: an intermediate QP along the homotopy path was solved.\n" );
			break;

		case QPS_SOLVED:
			myPrintf( "Status of QP object: solution of the actual QP was found.\n" );
			break;
	}

	switch ( options.printLevel )
	{
		case PL_LOW:
					myPrintf( "Print level of QP object is low, i.e. only error are printed.\n" );
			break;

		case PL_MEDIUM:
			myPrintf( "Print level of QP object is medium, i.e. error and warnings are printed.\n" );
			break;

		case PL_HIGH:
			myPrintf( "Print level of QP object is high, i.e. all available output is printed.\n" );
			break;

		default:
			break;
	}

	myPrintf( "\n" );
	#endif
	#endif

	return SUCCESSFUL_RETURN;
}


returnValue QProblemB::printOptions( ) const
{
	return options.print( );
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue QProblemB::clear( )
{
	if ( ( freeHessian == BT_TRUE ) && ( H != 0 ) )
	{
		delete H;
		H = 0;
	}

	if ( g != 0 )
	{
		delete[] g;
		g = 0;
	}

	if ( lb != 0 )
	{
		delete[] lb;
		lb = 0;
	}

	if ( ub != 0 )
	{
		delete[] ub;
		ub = 0;
	}

	if ( R != 0 )
	{
		delete[] R;
		R = 0;
	}

	if ( x != 0 )
	{
		delete[] x;
		x = 0;
	}

	if ( y != 0 )
	{
		delete[] y;
		y = 0;
	}

	if ( delta_xFR_TMP != 0 )
	{
		delete[] delta_xFR_TMP;
		delta_xFR_TMP = 0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue QProblemB::copy(	const QProblemB& rhs
								)
{
	int _nV = rhs.getNV( );

	bounds = rhs.bounds;

	freeHessian = rhs.freeHessian;
	
	if ( freeHessian == BT_TRUE )
		H = dynamic_cast<SymmetricMatrix *>(rhs.H->duplicate());
	else
		H = rhs.H;

	if ( rhs.g != 0 )
	{
		g = new real_t[_nV];
		setG( rhs.g );
	}
	else
		g = 0;

	if ( rhs.lb != 0 )
	{
		lb = new real_t[_nV];
		setLB( rhs.lb );
	}
	else
		lb = 0;

	if ( rhs.ub != 0 )
	{
		ub = new real_t[_nV];
		setUB( rhs.ub );
	}
	else
		ub = 0;

	if ( rhs.R != 0 )
	{
		R = new real_t[_nV*_nV];
		memcpy( R,rhs.R,_nV*_nV*sizeof(real_t) );
	}
	else
		R = 0;
	
	haveCholesky = rhs.haveCholesky;

	if ( rhs.x != 0 )
	{
		x = new real_t[_nV];
		memcpy( x,rhs.x,_nV*sizeof(real_t) );
	}
	else
		x = 0;

	if ( rhs.y != 0 )
	{
		y = new real_t[_nV];
		memcpy( y,rhs.y,_nV*sizeof(real_t) );
	}
	else
		y = 0;

	tau = rhs.tau;

	hessianType = rhs.hessianType;
	isRegularised = rhs.isRegularised;

	infeasible = rhs.infeasible;
	unbounded = rhs.unbounded;

	status = rhs.status;

	count = rhs.count;

	ramp0 = rhs.ramp0;
	ramp1 = rhs.ramp1;

	delta_xFR_TMP = new real_t[_nV];	/* nFR */

	options = rhs.options;
	setPrintLevel( options.printLevel );

	flipper = rhs.flipper;

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e H e s s i a n T y p e
 */
returnValue QProblemB::determineHessianType( )
{
	int i;
	int nV = getNV( );

	/* if Hessian type has been set by user, do NOT change it! */
	if ( hessianType != HST_UNKNOWN )
		return SUCCESSFUL_RETURN;

	/* if Hessian has not been allocated, assume it to be all zeros! */
	if ( H == 0 )
	{
		hessianType = HST_ZERO;
		return SUCCESSFUL_RETURN;
	}


	/* 1) If Hessian has outer-diagonal elements,
	 *    Hessian is assumed to be positive definite. */
	hessianType = HST_POSDEF;
	if (H->isDiag() == BT_FALSE)
		return SUCCESSFUL_RETURN;

	/* 2) Otherwise it is diagonal and test for identity or zero matrix is performed. */
	/* hessianType = HST_DIAGONAL; */

	BooleanType isIdentity = BT_TRUE;
	BooleanType isZero = BT_TRUE;

	for ( i=0; i<nV; ++i )
	{
		if ( options.enableFlippingBounds == BT_FALSE )
			if ( H->diag(i) < -ZERO ) 
				return THROWERROR( RET_HESSIAN_INDEFINITE );

		if ( getAbs( H->diag(i) - 1.0 ) > EPS )
			isIdentity = BT_FALSE;

		if ( getAbs( H->diag(i) ) > EPS )
			isZero = BT_FALSE;
	}

	if ( isIdentity == BT_TRUE )
		hessianType = HST_IDENTITY;

	if ( isZero == BT_TRUE )
		hessianType = HST_ZERO;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblemB::setupSubjectToType( )
{
	return setupSubjectToType( lb,ub );
}


/*
 *	s e t u p S u b j e c t T o T y p e
 */
returnValue QProblemB::setupSubjectToType( const real_t* const lb_new, const real_t* const ub_new )
{
	int i;
	int nV = getNV( );


	/* 1) Check if lower bounds are present. */
	bounds.setNoLower( BT_TRUE );
	if ( lb_new != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( lb_new[i] > -INFTY )
			{
				bounds.setNoLower( BT_FALSE );
				break;
			}
		}
	}

	/* 2) Check if upper bounds are present. */
	bounds.setNoUpper( BT_TRUE );
	if ( ub_new != 0 )
	{
		for( i=0; i<nV; ++i )
		{
			if ( ub_new[i] < INFTY )
			{
				bounds.setNoUpper( BT_FALSE );
				break;
			}
		}
	}

	/* 3) Determine implicitly fixed and unbounded variables. */
	if ( ( lb_new != 0 ) && ( ub_new != 0 ) )
	{
		for( i=0; i<nV; ++i )
		{
			if ( ( lb_new[i] < -INFTY + options.boundTolerance ) && ( ub_new[i] > INFTY - options.boundTolerance )
					&& (options.enableFarBounds == BT_FALSE))
			{
				bounds.setType( i,ST_UNBOUNDED );
			}
			else
			{
				if ( options.enableEqualities && lb_new[i] > ub_new[i] - options.boundTolerance )
					bounds.setType( i,ST_EQUALITY );
				else
					bounds.setType( i,ST_BOUNDED );
			}
		}
	}
	else
	{
		if ( ( lb_new == 0 ) && ( ub_new == 0 ) )
		{
			for( i=0; i<nV; ++i )
				bounds.setType( i,ST_UNBOUNDED );
		}
		else
		{
			for( i=0; i<nV; ++i )
				bounds.setType( i,ST_BOUNDED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p C h o l e s k y D e c o m p o s i t i o n
 */
returnValue QProblemB::setupCholeskyDecomposition( )
{
	int i, j;
	int nV  = getNV( );
	int nFR = getNFR( );
	
	/* 1) Initialises R with all zeros. */
	for( i=0; i<nV*nV; ++i )
		R[i] = 0.0;

	/* 2) Calculate Cholesky decomposition of H (projected to free variables). */
	if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_IDENTITY ) )
	{
		if ( hessianType == HST_ZERO )
		{
			/* if Hessian is zero matrix, it is assumed that it has been
			 * regularised and thus its Cholesky factor is the identity
			 * matrix scaled by sqrt(eps). */
			if ( usingRegularisation( ) == BT_TRUE )
			{
				for( i=0; i<nV; ++i )
					RR(i,i) = sqrt( options.epsRegularisation );
			}
			else
				return THROWERROR( RET_CHOLESKY_OF_ZERO_HESSIAN );
		}
		else
		{
			/* if Hessian is identity, so is its Cholesky factor. */
			for( i=0; i<nV; ++i )
				RR(i,i) = 1.0;
		}
	}
	else
	{
		if ( nFR > 0 )
		{
			int* FR_idx;
			bounds.getFree( )->getNumberArray( &FR_idx );

			/* get H */
			for ( j=0; j < nFR; ++j )
				H->getCol (FR_idx[j], bounds.getFree (), 1.0, &R[j*nV]);

			/* R'*R = H */
			long info = 0;
			unsigned long _nFR = nFR, _nV = nV;

			POTRF ("U", &_nFR, R, &_nV, &info);

			/* <0 = invalid call, =0 ok, >0 not spd */
			if (info > 0) {
				hessianType = HST_SEMIDEF;
				return RET_HESSIAN_NOT_SPD;
			}

			/* zero first subdiagonal to make givens updates work */
			for (i=0;i<nFR-1;++i)
				RR(i+1,i) = 0.0;

		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	o b t a i n A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB::obtainAuxiliaryWorkingSet(	const real_t* const xOpt, const real_t* const yOpt,
													const Bounds* const guessedBounds, Bounds* auxiliaryBounds
													) const
{
	int i = 0;
	int nV = getNV( );


	/* 1) Ensure that desiredBounds is allocated (and different from guessedBounds). */
	if ( ( auxiliaryBounds == 0 ) || ( auxiliaryBounds == guessedBounds ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 2) Setup working set for auxiliary initial QP. */
	if ( guessedBounds != 0 )
	{
		/* If an initial working set is specific, use it!
		 * Moreover, add all implictly fixed variables if specified. */
		for( i=0; i<nV; ++i )
		{
			#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
			if ( bounds.getType( i ) == ST_EQUALITY )
			{
				if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
			else
			#endif
			{
				if ( auxiliaryBounds->setupBound( i,guessedBounds->getStatus( i ) ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
			}
		}
	}
	else	/* No initial working set specified. */
	{
		if ( ( xOpt != 0 ) && ( yOpt == 0 ) )
		{
			/* Obtain initial working set by "clipping". */
			for( i=0; i<nV; ++i )
			{
				if ( xOpt[i] <= lb[i] + options.boundTolerance )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( xOpt[i] >= ub[i] - options.boundTolerance )
				{
					if ( auxiliaryBounds->setupBound( i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all implictly fixed variables if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( auxiliaryBounds->setupBound( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		if ( ( xOpt == 0 ) && ( yOpt != 0 ) )
		{
			/* Obtain initial working set in accordance to sign of dual solution vector. */
			for( i=0; i<nV; ++i )
			{
				if ( yOpt[i] > EPS )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				if ( yOpt[i] < -EPS )
				{
					if ( auxiliaryBounds->setupBound( i,ST_UPPER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
					continue;
				}

				/* Moreover, add all implictly fixed variables if specified. */
				#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
				else
				#endif
				{
					if ( auxiliaryBounds->setupBound( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
						return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
				}
			}
		}

		/* If xOpt and yOpt are null pointer and no initial working is specified,
		 * start with empty working set (or implicitly fixed bounds only)
		 * for auxiliary QP. */
		if ( ( xOpt == 0 ) && ( yOpt == 0 ) )
		{
			for( i=0; i<nV; ++i )
			{
				switch( bounds.getType( i ) )
				{
					case ST_UNBOUNDED:
						if ( auxiliaryBounds->setupBound( i,ST_INACTIVE ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;

					/* Only add all implictly fixed variables if specified. */
					#ifdef __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__
					case ST_EQUALITY:
						if ( auxiliaryBounds->setupBound( i,ST_LOWER ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;
					#endif

					default:
						if ( auxiliaryBounds->setupBound( i,options.initialStatusBounds ) != SUCCESSFUL_RETURN )
							return THROWERROR( RET_OBTAINING_WORKINGSET_FAILED );
						break;
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	b a c k s o l v e R
 */
returnValue QProblemB::backsolveR(	const real_t* const b, BooleanType transposed,
									real_t* const a
									) const
{
	/* Call standard backsolve procedure (i.e. removingBound == BT_FALSE). */
	return backsolveR( b,transposed,BT_FALSE,a );
}


/*
 *	b a c k s o l v e R
 */
returnValue QProblemB::backsolveR(	const real_t* const b, BooleanType transposed,
									BooleanType removingBound,
									real_t* const a
									) const
{
	int i, j;
	int nV = getNV( );
	int nR = getNZ( );

	real_t sum;

	/* if backsolve is called while removing a bound, reduce nZ by one. */
	if ( removingBound == BT_TRUE )
		--nR;

	/* nothing to do */
	if ( nR <= 0 )
		return SUCCESSFUL_RETURN;


	/* Solve Ra = b, where R might be transposed. */
	if ( transposed == BT_FALSE )
	{
		/* solve Ra = b */
		for( i=(nR-1); i>=0; --i )
		{
			sum = b[i];
			for( j=(i+1); j<nR; ++j )
				sum -= RR(i,j) * a[j];

			if ( getAbs( RR(i,i) ) >= ZERO*getAbs( sum ) )
				a[i] = sum / RR(i,i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}
	else
	{
		/* solve R^T*a = b */
		for( i=0; i<nR; ++i )
		{
			sum = b[i];
			for( j=0; j<i; ++j )
				sum -= RR(j,i) * a[j];

			if ( getAbs( RR(i,i) ) >= ZERO*getAbs( sum ) )
				a[i] = sum / RR(i,i);
			else
				return THROWERROR( RET_DIV_BY_ZERO );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e D a t a S h i f t
 */
returnValue QProblemB::determineDataShift(	const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new,
											real_t* const delta_g, real_t* const delta_lb, real_t* const delta_ub,
											BooleanType& Delta_bB_isZero
											)
{
	int i, ii;
	int nV  = getNV( );
	int nFX = getNFX( );

	int* FX_idx;
	bounds.getFixed( )->getNumberArray( &FX_idx );


	/* 1) Calculate shift directions. */
	for( i=0; i<nV; ++i )
		delta_g[i]  = g_new[i]  - g[i];

	if ( lb_new != 0 )
	{
		for( i=0; i<nV; ++i )
			delta_lb[i] = lb_new[i] - lb[i];
	}
	else
	{
		/* if no lower bounds exist, assume the new lower bounds to be -infinity */
		for( i=0; i<nV; ++i )
			delta_lb[i] = -INFTY - lb[i];
	}

	if ( ub_new != 0 )
	{
		for( i=0; i<nV; ++i )
			delta_ub[i] = ub_new[i] - ub[i];
	}
	else
	{
		/* if no upper bounds exist, assume the new upper bounds to be infinity */
		for( i=0; i<nV; ++i )
			delta_ub[i] = INFTY - ub[i];
	}

	/* 2) Determine if active bounds are to be shifted. */
	Delta_bB_isZero = BT_TRUE;

	for ( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];

		if ( ( getAbs( delta_lb[ii] ) > EPS ) || ( getAbs( delta_ub[ii] ) > EPS ) )
		{
			Delta_bB_isZero = BT_FALSE;
			break;
		}
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p Q P d a t a
 */
returnValue QProblemB::setupQPdata(	SymmetricMatrix *_H, const real_t* const _g,
									const real_t* const _lb, const real_t* const _ub
									)
{
	int i;
	int nV = getNV( );

	/* 1) Setup Hessian matrix. */
	setH( _H );

	/* 2) Setup gradient vector. */
	if ( _g == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );
	else
		setG( _g );

	/* 3) Setup lower bounds vector. */
	if ( _lb != 0 )
		setLB( _lb );
	else
		/* if no lower bounds are specified, set them to -infinity */
		for( i=0; i<nV; ++i )
			lb[i] = -INFTY;

	/* 4) Setup upper bounds vector. */
	if ( _ub != 0 )
		setUB( _ub );
	else
		/* if no upper bounds are specified, set them to infinity */
		for( i=0; i<nV; ++i )
			ub[i] = INFTY;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a
 */
returnValue QProblemB::setupQPdata(	const real_t* const _H, const real_t* const _g,
									const real_t* const _lb, const real_t* const _ub
									)
{
	int i;
	int nV = getNV( );

	/* 1) Setup Hessian matrix. */
	setH( _H );

	/* 2) Setup gradient vector. */
	if ( _g == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );
	else
		setG( _g );

	/* 3) Setup lower bounds vector. */
	if ( _lb != 0 )
	{
		setLB( _lb );
	}
	else
	{
		/* if no lower bounds are specified, set them to -infinity */
		for( i=0; i<nV; ++i )
			lb[i] = -INFTY;
	}

	/* 4) Setup upper bounds vector. */
	if ( _ub != 0 )
	{
		setUB( _ub );
	}
	else
	{
		/* if no upper bounds are specified, set them to infinity */
		for( i=0; i<nV; ++i )
			ub[i] = INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p Q P d a t a F r o m F i l e
 */
returnValue QProblemB::setupQPdataFromFile(	const char* const H_file, const char* const g_file,
											const char* const lb_file, const char* const ub_file
											)
{
	int i;
	int nV = getNV( );

	returnValue returnvalue;


	/* 1) Load Hessian matrix from file. */
	if ( H_file != 0 )
	{
		real_t* _H = new real_t[nV * nV];
		returnvalue = readFromFile( _H, nV,nV, H_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] _H;
			return THROWERROR( returnvalue );
		}
		setH( _H );
		H->doFreeMemory( );
	}
	else
	{
		real_t* _H = 0;
		setH( _H );
	}

	/* 2) Load gradient vector from file. */
	if ( g_file == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	returnvalue = readFromFile( g, nV, g_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
		return THROWERROR( returnvalue );

	/* 3) Load lower bounds vector from file. */
	if ( lb_file != 0 )
	{
		returnvalue = readFromFile( lb, nV, lb_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* if no lower bounds are specified, set them to -infinity */
		for( i=0; i<nV; ++i )
			lb[i] = -INFTY;
	}

	/* 4) Load upper bounds vector from file. */
	if ( ub_file != 0 )
	{
		returnvalue = readFromFile( ub, nV, ub_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* if no upper bounds are specified, set them to infinity */
		for( i=0; i<nV; ++i )
			ub[i] = INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	l o a d Q P v e c t o r s F r o m F i l e
 */
returnValue QProblemB::loadQPvectorsFromFile(	const char* const g_file, const char* const lb_file, const char* const ub_file,
												real_t* const g_new, real_t* const lb_new, real_t* const ub_new
												) const
{
	int nV = getNV( );

	returnValue returnvalue;


	/* 1) Load gradient vector from file. */
	if ( ( g_file != 0 ) && ( g_new != 0 ) )
	{
		returnvalue = readFromFile( g_new, nV, g_file );
		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		/* At least gradient vector needs to be specified! */
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	/* 2) Load lower bounds vector from file. */
	if ( lb_file != 0 )
	{
		if ( lb_new != 0 )
		{
			returnvalue = readFromFile( lb_new, nV, lb_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* If filename is given, storage must be provided! */
			return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	/* 3) Load upper bounds vector from file. */
	if ( ub_file != 0 )
	{
		if ( ub_new != 0 )
		{
			returnvalue = readFromFile( ub_new, nV, ub_file );
			if ( returnvalue != SUCCESSFUL_RETURN )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* If filename is given, storage must be provided! */
			return THROWERROR( RET_INVALID_ARGUMENTS );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t I n f e a s i b i l i t y F l a g
 */
returnValue QProblemB::setInfeasibilityFlag(	returnValue returnvalue
												)
{
	infeasible = BT_TRUE;

	if ( options.enableFarBounds == BT_FALSE )
		THROWERROR( returnvalue );

	return returnvalue;
}


/*
 *	i s C P U t i m e L i m i t E x c e e d e d
 */
BooleanType QProblemB::isCPUtimeLimitExceeded(	const real_t* const cputime,
												real_t starttime,
												int nWSR
												) const
{
	/* Always perform next QP iteration if no CPU time limit is given. */
	if ( cputime == 0 )
		return BT_FALSE;

	/* Always perform first QP iteration. */
	if ( nWSR <= 0 )
		return BT_FALSE;

	real_t elapsedTime = getCPUtime( ) - starttime;
	real_t timePerIteration = elapsedTime / ((real_t) nWSR);

	/* Determine if next QP iteration exceed CPU time limit
	 * considering the (current) average CPU time per iteration. */
	if ( ( elapsedTime + timePerIteration*1.25 ) <= ( *cputime ) )
		return BT_FALSE;
	else
		return BT_TRUE;
}


/*
 *	r e g u l a r i s e H e s s i a n
 */
returnValue QProblemB::regulariseHessian( )
{
	/* Do nothing if Hessian regularisation is disbaled! */
	if ( options.enableRegularisation == BT_FALSE )
		return SUCCESSFUL_RETURN;

	/* Regularisation of identity Hessian not possible. */
	if ( hessianType == HST_IDENTITY )
		return THROWERROR( RET_CANNOT_REGULARISE_IDENTITY );

	/* Determine regularisation parameter. */
	if ( usingRegularisation( ) == BT_TRUE )
		return THROWERROR( RET_HESSIAN_ALREADY_REGULARISED );
	else
	{
		/* Regularisation of zero Hessian is done implicitly. */
		if ( hessianType != HST_ZERO )
			if ( H->addToDiag( H->getNorm() * options.epsRegularisation ) == RET_NO_DIAGONAL_AVAILABLE )
				return THROWERROR( RET_CANNOT_REGULARISE_SPARSE );

		isRegularised = BT_TRUE;
		THROWINFO( RET_USING_REGULARISATION );
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	p e r f o r m R a t i o T e s t
 */
returnValue QProblemB::performRatioTest(	const int nIdx,
											const int* const idxList,
											const SubjectTo* const subjectTo,
											const real_t* const num,
											const real_t* const den,
											real_t epsNum,
											real_t epsDen,
											real_t& t,
											int& BC_idx
											) const
{
	int i, ii;

	BC_idx = -1;

	for( i=0; i<nIdx; ++i )
	{
		ii = idxList[i];

		if ( subjectTo->getType( ii ) != ST_EQUALITY )
		{
			if ( ( subjectTo->getStatus( ii ) == ST_LOWER ) || ( subjectTo->getStatus( ii ) == ST_INACTIVE ) )
			{
				if ( isBlocking( num[i],den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					t = num[i] / den[i];
					BC_idx = ii;
				}
			}
			else
			if ( subjectTo->getStatus( ii ) == ST_UPPER )
			{
				if ( isBlocking( -num[i],-den[i],epsNum,epsDen,t ) == BT_TRUE )
				{
					t = num[i] / den[i];
					BC_idx = ii;
				}
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 * r e l a t i v e H o m o t o p y L e n g t h
 */
real_t QProblemB::relativeHomotopyLength(const real_t* const g_new, const real_t* const lb_new, const real_t* const ub_new)
{
	int nV = getNV( ), i;
	real_t d, s, len = 0.0;

	/* gradient */
	for (i = 0; i < nV; i++)
	{
		s = getAbs(g_new[i]);
		if (s < 1.0) s = 1.0;
		d = getAbs(g_new[i] - g[i]) / s;
		if (d > len) len = d;
	}

	/* lower bounds */
	for (i = 0; i < nV && lb_new; i++)
	{
		s = getAbs(lb_new[i]);
		if (s < 1.0) s = 1.0;
		d = getAbs(lb_new[i] - lb[i]) / s;
		if (d > len) len = d;
	}

	/* upper bounds */
	for (i = 0; i < nV && ub_new; i++)
	{
		s = getAbs(ub_new[i]);
		if (s < 1.0) s = 1.0;
		d = getAbs(ub_new[i] - ub[i]) / s;
		if (d > len) len = d;
	}

	return len;
}


/*
 * p e r f o r m R a m p i n g
 */
returnValue QProblemB::performRamping( )
{
	int nV = getNV( ), bstat, i;
	real_t t, rampVal;

	/* ramp inactive bounds and active dual variables */
	for (i = 0; i < nV; i++)
	{
		switch (bounds.getType(i))
		{
			case ST_EQUALITY: lb[i] = x[i]; ub[i] = x[i]; continue; /* reestablish exact feasibility */
			case ST_UNBOUNDED: continue;
			case ST_DISABLED: continue;
			default: break;
		}

		t = static_cast<real_t>((i + rampOffset) % nV) / static_cast<real_t>(nV-1);
		rampVal = (1.0-t) * ramp0 + t * ramp1;
		bstat = bounds.getStatus(i);
		if (bstat != ST_LOWER) { lb[i] = x[i] - rampVal; }
		if (bstat != ST_UPPER) { ub[i] = x[i] + rampVal; }
		if (bstat == ST_LOWER) { lb[i] = x[i]; y[i] = +rampVal; }
		if (bstat == ST_UPPER) { ub[i] = x[i]; y[i] = -rampVal; }
		if (bstat == ST_INACTIVE) y[i] = 0.0; /* reestablish exact complementarity */
	}

	/* reestablish exact stationarity */
	setupAuxiliaryQPgradient( );

	/* advance ramp offset to avoid Ramping cycles */
	rampOffset++;

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R I V A T E                                                            *
 *****************************************************************************/

/*
 *	s o l v e I n i t i a l Q P
 */
returnValue QProblemB::solveInitialQP(	const real_t* const xOpt, const real_t* const yOpt,
										const Bounds* const guessedBounds,
										int& nWSR, real_t* const cputime
										)
{
	int i;
	int nV = getNV( );


	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = getCPUtime( );


	status = QPS_NOTINITIALISED;

	/* I) ANALYSE QP DATA: */
	/* 1) Check if Hessian happens to be the identity matrix. */
	if ( determineHessianType( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

/*	char myPrintfString[40];
	snprintf( myPrintfString,40,"\ntype: %d\n",hessianType );
	myPrintf( myPrintfString );*/
	
	/* 2) Setup type of bounds (i.e. unbounded, implicitly fixed etc.). */
	if ( setupSubjectToType( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	status = QPS_PREPARINGAUXILIARYQP;


	/* II) SETUP AUXILIARY QP WITH GIVEN OPTIMAL SOLUTION: */
	/* 1) Setup bounds data structure. */
	if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 2) Setup optimal primal/dual solution for auxiliary QP. */
	if ( setupAuxiliaryQPsolution( xOpt,yOpt ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 3) Obtain linear independent working set for auxiliary QP. */
	Bounds auxiliaryBounds( nV );
	if ( obtainAuxiliaryWorkingSet( xOpt,yOpt,guessedBounds, &auxiliaryBounds ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* 4) Setup working set of auxiliary QP and setup cholesky decomposition. */
	/* a) Working set of auxiliary QP. */
	if ( setupAuxiliaryWorkingSet( &auxiliaryBounds,BT_TRUE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_INIT_FAILED );

	/* b) Regularise Hessian if necessary. */
	if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_SEMIDEF ) )
	{
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_INIT_FAILED_REGULARISATION );
	}

	haveCholesky = BT_FALSE;

	/* c) Cholesky decomposition. */
	returnValue returnvalueCholesky = setupCholeskyDecomposition( );

	/* If Hessian is not positive definite, regularise and try again. */
	if ( returnvalueCholesky == RET_HESSIAN_NOT_SPD )
	{
		if ( regulariseHessian( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_INIT_FAILED );

		if ( setupCholeskyDecomposition( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_INIT_FAILED_CHOLESKY );
	}

	if ( returnvalueCholesky == RET_INDEXLIST_CORRUPTED )
		return THROWERROR( RET_INIT_FAILED_CHOLESKY );

	haveCholesky = BT_TRUE;

	/* 5) Store original QP formulation... */
	real_t* g_original = new real_t[nV];
	real_t* lb_original = new real_t[nV];
	real_t* ub_original = new real_t[nV];

	for( i=0; i<nV; ++i )
	{
		g_original[i]  = g[i];
		lb_original[i] = lb[i];
		ub_original[i] = ub[i];
	}

	/* ... and setup QP data of an auxiliary QP having an optimal solution
	 * as specified by the user (or xOpt = yOpt = 0, by default). */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
	{
		delete[] ub_original; delete[] lb_original; delete[] g_original;
		return THROWERROR( RET_INIT_FAILED );
	}

	if ( setupAuxiliaryQPbounds( BT_TRUE ) != SUCCESSFUL_RETURN )
	{
 		delete[] ub_original; delete[] lb_original; delete[] g_original;
		return THROWERROR( RET_INIT_FAILED );
	}

	status = QPS_AUXILIARYQPSOLVED;


	/* III) SOLVE ACTUAL INITIAL QP: */

	/* Allow only remaining CPU time for usual hotstart. */
	if ( cputime != 0 )
		*cputime -= getCPUtime( ) - starttime;

	/* Use hotstart method to find the solution of the original initial QP,... */
	returnValue returnvalue = hotstart( g_original,lb_original,ub_original, nWSR,cputime );

	/* ... deallocate memory,... */
	delete[] ub_original; delete[] lb_original; delete[] g_original;

	/* ... check for infeasibility and unboundedness... */
	if ( isInfeasible( ) == BT_TRUE )
		return THROWERROR( RET_INIT_FAILED_INFEASIBILITY );

	if ( isUnbounded( ) == BT_TRUE )
		return THROWERROR( RET_INIT_FAILED_UNBOUNDEDNESS );

	/* ... and internal errors. */
	if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
		return THROWERROR( RET_INIT_FAILED_HOTSTART );


	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = getCPUtime( ) - starttime;

	THROWINFO( RET_INIT_SUCCESSFUL );

	return returnvalue;
}


/*
 *	s o l v e Q P
 */
returnValue QProblemB::solveQP(	const real_t* const g_new,
								const real_t* const lb_new, const real_t* const ub_new,
								int& nWSR, real_t* const cputime, int nWSRperformed
								)
{
	int iter;
	int nV  = getNV( );

	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )       ||
		 ( getStatus( ) == QPS_PREPARINGAUXILIARYQP ) ||
		 ( getStatus( ) == QPS_PERFORMINGHOMOTOPY )   )
	{
		return THROWERROR( RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED );
	}

	/* start runtime measurement */
	real_t starttime = 0.0;
	if ( cputime != 0 )
		starttime = getCPUtime( );


	/* I) PREPARATIONS */
	/* 1) Allocate delta vectors of gradient and bounds,
	 *    index arrays and step direction arrays. */
	int nFR, nFX;

	real_t* delta_xFR = new real_t[nV];
	real_t* delta_xFX = new real_t[nV];
	real_t* delta_yFX = new real_t[nV];

	real_t* delta_g  = new real_t[nV];
	real_t* delta_lb = new real_t[nV];
	real_t* delta_ub = new real_t[nV];

	returnValue returnvalue;
	BooleanType Delta_bB_isZero;

	int BC_idx;
	SubjectToStatus BC_status;

	real_t homotopyLength;

	char messageString[80];

	/* 2) Update type of bounds, e.g. a formerly implicitly fixed
	 *    variable might have become a normal one etc. */
	if ( setupSubjectToType( lb_new,ub_new ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_HOTSTART_FAILED );

	/* 3) Reset status flags. */
	infeasible = BT_FALSE;
	unbounded  = BT_FALSE;


	/* II) MAIN HOMOTOPY LOOP */
	for( iter=nWSRperformed; iter<nWSR; ++iter )
	{
		if ( isCPUtimeLimitExceeded( cputime,starttime,iter-nWSRperformed ) == BT_TRUE )
		{
			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			break;
		}

		status = QPS_PERFORMINGHOMOTOPY;

		#ifndef __XPCTARGET__
		snprintf( messageString,80,"%d ...",iter );
		getGlobalMessageHandler( )->throwInfo( RET_ITERATION_STARTED,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		#endif

		/* some more definitions */
		nFR = getNFR( );
		nFX = getNFX( );

		/* 2) Initialise shift direction of the gradient and the bounds. */
		returnvalue = determineDataShift(	g_new,lb_new,ub_new,
											delta_g,delta_lb,delta_ub,
											Delta_bB_isZero
											);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_SHIFT_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 3) Determination of step direction of X and Y. */
		returnvalue = determineStepDirection(	delta_g,delta_lb,delta_ub,
												Delta_bB_isZero,
												delta_xFX,delta_xFR,delta_yFX
												);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_STEPDIRECTION_DETERMINATION_FAILED );
			return returnvalue;
		}


		/* 4) Determination of step length TAU.
		 *    This step along the homotopy path is also taken (without changing working set). */
		returnvalue = performStep(	delta_g,delta_lb,delta_ub,
									delta_xFX,delta_xFR,delta_yFX,
									BC_idx,BC_status
									);
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			THROWERROR( RET_STEPLENGTH_DETERMINATION_FAILED );
			return returnvalue;
		}

		/* 5) Termination criterion. */
		homotopyLength = relativeHomotopyLength(g_new, lb_new, ub_new);
		if ( homotopyLength <= options.terminationTolerance )
		{
			status = QPS_SOLVED;

			THROWINFO( RET_OPTIMAL_SOLUTION_FOUND );

			if ( options.printLevel > PL_NONE )
				if ( printIteration( iter,BC_idx,BC_status ) != SUCCESSFUL_RETURN )
					THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass this as return value! */

			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			return SUCCESSFUL_RETURN;
		}


		/* 6) Change active set. */
		returnvalue = changeActiveSet( BC_idx,BC_status );

		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
			delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

			/* Assign number of working set recalculations and stop runtime measurement. */
			nWSR = iter;
			if ( cputime != 0 )
				*cputime = getCPUtime( ) - starttime;

			/* checks for infeasibility... */
			if ( infeasible == BT_TRUE )
			{
				status = QPS_HOMOTOPYQPSOLVED;
				return setInfeasibilityFlag( RET_HOTSTART_STOPPED_INFEASIBILITY );
			}

			/* ...unboundedness... */
			if ( unbounded == BT_TRUE ) /* not necessary since objective function convex! */
				return THROWERROR( RET_HOTSTART_STOPPED_UNBOUNDEDNESS );

			/* ... and throw unspecific error otherwise */
			THROWERROR( RET_HOMOTOPY_STEP_FAILED );
			return returnvalue;
		}

		/* 7) Perform Ramping Strategy on zero homotopy step or drift correction (if desired). */
		 if ( ( tau <= EPS ) && ( options.enableRamping == BT_TRUE ) )
			performRamping( );
		else
		if ( (options.enableDriftCorrection > 0)
		  && ((iter+1) % options.enableDriftCorrection == 0) )
			performDriftCorrection( );  /* always returns SUCCESSFUL_RETURN */

		/* 8) Output information of successful QP iteration. */
		status = QPS_HOMOTOPYQPSOLVED;

		if ( options.printLevel > PL_NONE )
			if ( printIteration( iter,BC_idx,BC_status ) != SUCCESSFUL_RETURN )
				THROWERROR( RET_PRINT_ITERATION_FAILED ); /* do not pass this as return value! */
	}

	delete[] delta_yFX; delete[] delta_xFX; delete[] delta_xFR;
	delete[] delta_ub; delete[] delta_lb; delete[] delta_g;

	/* stop runtime measurement */
	if ( cputime != 0 )
		*cputime = getCPUtime( ) - starttime;


	/* if programm gets to here, output information that QP could not be solved
	 * within the given maximum numbers of working set changes */
	if ( options.printLevel == PL_HIGH )
	{
		#ifndef __XPCTARGET__
		snprintf( messageString,80,"(nWSR = %d)",iter );
		return getGlobalMessageHandler( )->throwWarning( RET_MAX_NWSR_REACHED,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		#else
		return RET_MAX_NWSR_REACHED;
		#endif
	}
	else
	{
		return RET_MAX_NWSR_REACHED;
	}
}


/*
 *	s o l v e R e g u l a r i s e d Q P
 */
returnValue QProblemB::solveRegularisedQP(	const real_t* const g_new,
											const real_t* const lb_new, const real_t* const ub_new,
											int& nWSR, real_t* const cputime, int nWSRperformed
											)
{
	int i, step;
	int nV = getNV( );


	/* Stop here if QP has not been regularised (i.e. normal QP solution). */
	if ( usingRegularisation( ) == BT_FALSE )
		return solveQP( g_new,lb_new,ub_new, nWSR,cputime,nWSRperformed );


	/* I) SOLVE USUAL REGULARISED QP */
	returnValue returnvalue;

	int nWSR_max   = nWSR;
	int nWSR_total = nWSRperformed;

	real_t cputime_total = 0.0;
	real_t cputime_cur   = 0.0;

	if ( cputime == 0 )
	{
		returnvalue = solveQP( g_new,lb_new,ub_new, nWSR,0,nWSRperformed );
	}
	else
	{
		cputime_cur = *cputime;
		returnvalue = solveQP( g_new,lb_new,ub_new, nWSR,&cputime_cur,nWSRperformed );
	}
	nWSR_total     = nWSR;
	cputime_total += cputime_cur;


	/* Only continue if QP solution has been successful. */
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		if ( cputime != 0 )
			*cputime = cputime_total;

		if ( returnvalue == RET_MAX_NWSR_REACHED )
			THROWWARNING( RET_NO_REGSTEP_NWSR );

		return returnvalue;
	}


	/* II) PERFORM SUCCESSIVE REGULARISATION STEPS */
	real_t* gMod = new real_t[nV];

	for( step=0; step<options.numRegularisationSteps; ++step )
	{
		/* 1) Modify gradient: gMod = g - eps*xOpt
		 *    (assuming regularisation matrix to be eps*Id). */
		for( i=0; i<nV; ++i )
			gMod[i] = g_new[i] - options.epsRegularisation*x[i];

		/* 2) Solve regularised QP with modified gradient allowing
		 *    only as many working set recalculations and CPU time
		 *    as have been left from previous QP solutions. */
		if ( cputime == 0 )
		{
			nWSR = nWSR_max;
			returnvalue = solveQP( gMod,lb_new,ub_new, nWSR,0,nWSR_total );
		}
		else
		{
			nWSR = nWSR_max;
			cputime_cur = *cputime - cputime_total;
			returnvalue = solveQP( gMod,lb_new,ub_new, nWSR,&cputime_cur,nWSR_total );
		}

		nWSR_total     = nWSR;
		cputime_total += cputime_cur;

		/* Only continue if QP solution has been successful. */
		if ( returnvalue != SUCCESSFUL_RETURN )
		{
			delete[] gMod;

			if ( cputime != 0 )
				*cputime = cputime_total;

			if ( returnvalue == RET_MAX_NWSR_REACHED )
				THROWWARNING( RET_FEWER_REGSTEPS_NWSR );

			return returnvalue;
		}
	}

	delete[] gMod;

	if ( cputime != 0 )
		*cputime = cputime_total;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB::setupAuxiliaryWorkingSet( 	const Bounds* const auxiliaryBounds,
													BooleanType setupAfresh
													)
{
	int i;
	int nV = getNV( );

	/* consistency checks */
	if ( auxiliaryBounds != 0 )
	{
		for( i=0; i<nV; ++i )
			if ( ( bounds.getStatus( i ) == ST_UNDEFINED ) || ( auxiliaryBounds->getStatus( i ) == ST_UNDEFINED ) )
				return THROWERROR( RET_UNKNOWN_BUG );
	}
	else
	{
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}


	/* I) SETUP CHOLESKY FLAG:
	 *    Cholesky decomposition shall only be updated if working set
	 *    shall be updated (i.e. NOT setup afresh!) */
	BooleanType updateCholesky;
	if ( setupAfresh == BT_TRUE )
		updateCholesky = BT_FALSE;
	else
		updateCholesky = BT_TRUE;


	/* II) REMOVE FORMERLY ACTIVE BOUNDS (IF NECESSARY): */
	if ( setupAfresh == BT_FALSE )
	{
		/* Remove all active bounds that shall be inactive AND
		*  all active bounds that are active at the wrong bound. */
		for( i=0; i<nV; ++i )
		{
			if ( ( bounds.getStatus( i ) == ST_LOWER ) && ( auxiliaryBounds->getStatus( i ) != ST_LOWER ) )
				if ( removeBound( i,updateCholesky ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );

			if ( ( bounds.getStatus( i ) == ST_UPPER ) && ( auxiliaryBounds->getStatus( i ) != ST_UPPER ) )
				if ( removeBound( i,updateCholesky ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}


	/* III) ADD NEWLY ACTIVE BOUNDS: */
	/*      Add all inactive bounds that shall be active AND
	 *      all formerly active bounds that have been active at the wrong bound. */
	for( i=0; i<nV; ++i )
	{
		if ( ( bounds.getStatus( i ) == ST_INACTIVE ) && ( auxiliaryBounds->getStatus( i ) != ST_INACTIVE ) )
		{
			if ( addBound( i,auxiliaryBounds->getStatus( i ),updateCholesky ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_WORKINGSET_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P s o l u t i o n
 */
returnValue QProblemB::setupAuxiliaryQPsolution(	const real_t* const xOpt, const real_t* const yOpt
													)
{
	int i;
	int nV = getNV( );


	/* Setup primal/dual solution vectors for auxiliary initial QP:
	 * if a null pointer is passed, a zero vector is assigned;
	 * old solution vector is kept if pointer to internal solution vector is passed. */
	if ( xOpt != 0 )
	{
		if ( xOpt != x )
			for( i=0; i<nV; ++i )
				x[i] = xOpt[i];
	}
	else
	{
		for( i=0; i<nV; ++i )
			x[i] = 0.0;
	}

	if ( yOpt != 0 )
	{
		if ( yOpt != y )
			for( i=0; i<nV; ++i )
				y[i] = yOpt[i];
	}
	else
	{
		for( i=0; i<nV; ++i )
			y[i] = 0.0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P g r a d i e n t
 */
returnValue QProblemB::setupAuxiliaryQPgradient( )
{
	int i;
	int nV = getNV( );

	/* Setup gradient vector: g = -H*x + y'*Id. */
	switch ( hessianType )
	{
		case HST_ZERO:
			if ( usingRegularisation( ) == BT_FALSE )
				for ( i=0; i<nV; ++i )
					g[i] = y[i];
			else
				for ( i=0; i<nV; ++i )
					g[i] = y[i] - options.epsRegularisation*x[i];
			break;

		case HST_IDENTITY:
			for ( i=0; i<nV; ++i )
				g[i] = y[i] - x[i];
			break;

		default:
			/* y'*Id */
			for ( i=0; i<nV; ++i )
				g[i] = y[i];

			/* -H*x */
			H->times(1, -1.0, x, nV, 1.0, g, nV);

			break;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P b o u n d s
 */
returnValue QProblemB::setupAuxiliaryQPbounds( BooleanType useRelaxation )
{
	int i;
	int nV = getNV( );


	/* Setup bound vectors. */
	for ( i=0; i<nV; ++i )
	{
		switch ( bounds.getStatus( i ) )
		{
			case ST_INACTIVE:
				if ( useRelaxation == BT_TRUE )
				{
					if ( bounds.getType( i ) == ST_EQUALITY )
					{
						lb[i] = x[i];
						ub[i] = x[i];
					}
					else
					{
						lb[i] = x[i] - options.boundRelaxation;
						ub[i] = x[i] + options.boundRelaxation;
					}
				}
				break;

			case ST_LOWER:
				lb[i] = x[i];
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					ub[i] = x[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						ub[i] = x[i] + options.boundRelaxation;
				}
				break;

			case ST_UPPER:
				ub[i] = x[i];
				if ( bounds.getType( i ) == ST_EQUALITY )
				{
					lb[i] = x[i];
				}
				else
				{
					if ( useRelaxation == BT_TRUE )
						lb[i] = x[i] - options.boundRelaxation;
				}
				break;

            case ST_DISABLED:
                break;
                
			default:
				return THROWERROR( RET_UNKNOWN_BUG );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A u x i l i a r y Q P
 */
returnValue QProblemB::setupAuxiliaryQP( const Bounds* const guessedBounds )
{
	int i;
	int nV = getNV( );

	/* nothing to do */
	if ( guessedBounds == &bounds )
		return SUCCESSFUL_RETURN;

	status = QPS_PREPARINGAUXILIARYQP;


	/* I) SETUP WORKING SET ... */
	if ( shallRefactorise( guessedBounds ) == BT_TRUE )
	{
		/* ... WITH REFACTORISATION: */
		/* 1) Reset bounds ... */
		bounds.init( nV );

		/*    ... and set them up afresh. */
		if ( setupSubjectToType( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		if ( bounds.setupAllFree( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 2) Setup guessed working set afresh. */
		if ( setupAuxiliaryWorkingSet( guessedBounds,BT_TRUE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

		/* 3) Calculate Cholesky decomposition. */
		if ( setupCholeskyDecomposition( ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}
	else
	{
		/* ... WITHOUT REFACTORISATION: */
		if ( setupAuxiliaryWorkingSet( guessedBounds,BT_FALSE ) != SUCCESSFUL_RETURN )
			THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );
	}


	/* II) SETUP AUXILIARY QP DATA: */
	/* 1) Ensure that dual variable is zero for fixed bounds. */
	for ( i=0; i<nV; ++i )
		if ( bounds.getStatus( i ) != ST_INACTIVE )
			y[i] = 0.0;

	/* 2) Setup gradient and bound vectors. */
	if ( setupAuxiliaryQPgradient( ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	if ( setupAuxiliaryQPbounds( BT_FALSE ) != SUCCESSFUL_RETURN )
		THROWERROR( RET_SETUP_AUXILIARYQP_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	d e t e r m i n e S t e p D i r e c t i o n
 */
returnValue QProblemB::determineStepDirection(	const real_t* const delta_g, const real_t* const delta_lb, const real_t* const delta_ub,
												BooleanType Delta_bB_isZero,
												real_t* const delta_xFX, real_t* const delta_xFR,
												real_t* const delta_yFX
												)
{
	int i, ii;
	int r;
	int nFR = getNFR( );
	int nFX = getNFX( );
	
	int* FR_idx;
	int* FX_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );


	/* This routine computes
	 * delta_xFX := delta_b
	 * delta_xFR := R \ R' \ -( delta_g + HMX*delta_xFX )
	 * delta_yFX := HMX'*delta_xFR + HFX*delta_xFX  { + eps*delta_xFX }
	 */

	/* I) DETERMINE delta_xFX := delta_{l|u}b */
	if ( Delta_bB_isZero == BT_FALSE )
	{
		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];

			if ( bounds.getStatus( ii ) == ST_LOWER )
				delta_xFX[i] = delta_lb[ii];
			else
				delta_xFX[i] = delta_ub[ii];
		}
	}
	else
	{
		for( i=0; i<nFX; ++i )
			delta_xFX[i] = 0.0;
	}


	/* delta_xFR_TMP holds the residual, initialized with right hand side
	 * delta_xFR holds the step that gets refined incrementally */
	for ( i=0; i<nFR; ++i )
	{
		ii = FR_idx[i];
		delta_xFR_TMP[i] = - delta_g[ii];
		delta_xFR[i] = 0.0;
	}


	/* Iterative refinement loop for delta_xFR */
	for ( r=0; r<=options.numRefinementSteps; ++r )
	{
		/* II) DETERMINE delta_xFR */
		if ( nFR > 0 )
		{
			/* Add - HMX*delta_xFX
			 * This is skipped if delta_b=0 or mixed part HM=0 (H=0 or H=Id) */
			if ( ( hessianType != HST_ZERO ) && ( hessianType != HST_IDENTITY ) && ( Delta_bB_isZero == BT_FALSE ) && ( r == 0 ) )
				H->times(bounds.getFree(), bounds.getFixed(), 1, -1.0, delta_xFX, nFX, 1.0, delta_xFR_TMP, nFR);

			/* Determine R' \ ( - HMX*delta_xFX - delta_gFR ) where R'R = HFR */
			if ( backsolveR( delta_xFR_TMP,BT_TRUE,delta_xFR_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );

			/* Determine HFR \ ( - HMX*delta_xFX - delta_gFR ) */
			if ( backsolveR( delta_xFR_TMP,BT_FALSE,delta_xFR_TMP ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_STEPDIRECTION_FAILED_CHOLESKY );
		}

		/* refine solution found for delta_xFR so far */
		for ( i=0; i<nFR; ++i )
			delta_xFR[i] += delta_xFR_TMP[i];

		if ( options.numRefinementSteps > 0 )
		{
			real_t rnrm = 0.0;
			/* compute new residual in delta_xFR_TMP:
			 * residual := - HFR*delta_xFR - HMX*delta_xFX - delta_gFR
			 * set to -delta_gFR */
			for ( i=0; i<nFR; ++i )
			{
				ii = FR_idx[i];
				delta_xFR_TMP[i] = -delta_g[ii];
			}
			/* add - HFR*delta_xFR */
			switch ( hessianType )
			{
				case HST_ZERO:
					break;

				case HST_IDENTITY:
					for ( i=0; i<nFR; ++i )
					{
						delta_xFR_TMP[i] -= delta_xFR[i];

						/* compute max norm */
						if (rnrm < getAbs (delta_xFR_TMP[i]))
							rnrm = getAbs (delta_xFR_TMP[i]);
					}
					break;

				default:
					H->times(bounds.getFree(), bounds.getFree(),  1, -1.0, delta_xFR, nFR, 1.0, delta_xFR_TMP, nFR);
					H->times(bounds.getFree(), bounds.getFixed(), 1, -1.0, delta_xFX, nFX, 1.0, delta_xFR_TMP, nFR);

					/* compute max norm */
					for ( i=0; i<nFR; ++i )
						if (rnrm < getAbs (delta_xFR_TMP[i]))
							rnrm = getAbs (delta_xFR_TMP[i]);

					break;
			}
			
			/* early termination of residual norm small enough */
			if ( rnrm < options.epsIterRef )
				break;
		}

	} /* end of refinement loop for delta_xFR */

	/* III) DETERMINE delta_yFX */
	if ( nFX > 0 )
	{
		if ( ( hessianType == HST_ZERO ) || ( hessianType == HST_IDENTITY ) )
		{
			for( i=0; i<nFX; ++i )
			{
				/* set to delta_g */
				ii = FX_idx[i];
				delta_yFX[i] = delta_g[ii];

				/* add HFX*delta_xFX = {0|I}*delta_xFX */
				if ( hessianType == HST_ZERO )
				{
					if ( usingRegularisation( ) == BT_TRUE )
						delta_yFX[i] += options.epsRegularisation*delta_xFX[i];
				}
				else
					delta_yFX[i] += delta_xFX[i];
			}
		}
		else
		{
			for( i=0; i<nFX; ++i )
			{
				/* set to delta_g */
				ii = FX_idx[i];
				delta_yFX[i] = delta_g[ii];
			}
			H->times(bounds.getFixed(), bounds.getFree(), 1, 1.0, delta_xFR, nFR, 1.0, delta_yFX, nFX);
			if (Delta_bB_isZero == BT_FALSE)
				H->times(bounds.getFixed(), bounds.getFixed(), 1, 1.0, delta_xFX, nFX, 1.0, delta_yFX, nFX);
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p e r f o r m S t e p
 */
returnValue QProblemB::performStep(	const real_t* const delta_g,
									const real_t* const delta_lb, const real_t* const delta_ub,
									const real_t* const delta_xFX,
									const real_t* const delta_xFR,
									const real_t* const delta_yFX,
									int& BC_idx, SubjectToStatus& BC_status
									)
{
	int i, ii;
	int nV = getNV( );
	int nFR = getNFR( );
	int nFX = getNFX( );

	int* FR_idx;
	int* FX_idx;

	bounds.getFree( )->getNumberArray( &FR_idx );
	bounds.getFixed( )->getNumberArray( &FX_idx );

	tau = 1.0;
	BC_idx = -1;
	BC_status = ST_UNDEFINED;

	int BC_idx_tmp = -1;

	real_t* num = new real_t[nV];
	real_t* den = new real_t[nV];


	/* I) DETERMINE MAXIMUM DUAL STEPLENGTH, i.e. ensure that
	 *    active dual bounds remain valid (ignoring implicitly fixed variables): */
	for( i=0; i<nFX; ++i )
	{
		ii = FX_idx[i];
		num[i] = y[ii];
		den[i] = -delta_yFX[i];
	}

	performRatioTest( nFX,FX_idx,&bounds, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

	if ( BC_idx_tmp >= 0 )
	{
		BC_idx = BC_idx_tmp;
		BC_status = ST_INACTIVE;
	}


	/* II) DETERMINE MAXIMUM PRIMAL STEPLENGTH, i.e. ensure that
	 *     inactive bounds remain valid (ignoring unbounded variables). */
	/* 1) Inactive lower bounds. */
	if ( bounds.hasNoLower( ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = getMax( x[ii] - lb[ii],0.0 );
			den[i] = delta_lb[ii] - delta_xFR[i];
		}

		performRatioTest( nFR,FR_idx,&bounds, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			BC_idx = BC_idx_tmp;
			BC_status = ST_LOWER;
		}
	}

	/* 2) Inactive upper bounds. */
	if ( bounds.hasNoUpper( ) == BT_FALSE )
	{
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			num[i] = getMax( ub[ii] - x[ii],0.0 );
			den[i] = delta_xFR[i] - delta_ub[ii];
		}

		performRatioTest( nFR,FR_idx,&bounds, num,den, options.epsNum,options.epsDen, tau,BC_idx_tmp );

		if ( BC_idx_tmp >= 0 )
		{
			BC_idx = BC_idx_tmp;
			BC_status = ST_UPPER;
		}
	}

	delete[] den;
	delete[] num;


	#ifndef __XPCTARGET__
	char messageString[80];

	if ( BC_status == ST_UNDEFINED )
		snprintf( messageString,80,"Stepsize is %.10e!",tau );
	else
		snprintf( messageString,80,"Stepsize is %.10e! (BC_idx = %d, BC_status = %d)",tau,BC_idx,BC_status );

	getGlobalMessageHandler( )->throwInfo( RET_STEPSIZE_NONPOSITIVE,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
	#endif


	/* III) PERFORM STEP ALONG HOMOTOPY PATH */
	if ( tau > ZERO )
	{
		/* 1) Perform step in primal und dual space. */
		for( i=0; i<nFR; ++i )
		{
			ii = FR_idx[i];
			x[ii] += tau*delta_xFR[i];
		}

		for( i=0; i<nFX; ++i )
		{
			ii = FX_idx[i];
			x[ii] += tau*delta_xFX[i];
			y[ii] += tau*delta_yFX[i];
		}

		/* 2) Shift QP data. */
		for( i=0; i<nV; ++i )
		{
			g[i]  += tau*delta_g[i];
			lb[i] += tau*delta_lb[i];
			ub[i] += tau*delta_ub[i];
		}
	}
	else
	{
		/* print a warning if stepsize is zero */
		#ifndef __XPCTARGET__
		snprintf( messageString,80,"Stepsize is %.6e",tau );
		getGlobalMessageHandler( )->throwWarning( RET_STEPSIZE,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
		#endif
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	c h a n g e A c t i v e S e t
 */
returnValue QProblemB::changeActiveSet( int BC_idx, SubjectToStatus BC_status )
{
	char messageString[80];

	/* IV) UPDATE ACTIVE SET */
	switch ( BC_status )
	{
		/* Optimal solution found as no working set change detected. */
		case ST_UNDEFINED:
			return RET_OPTIMAL_SOLUTION_FOUND;


		/* Remove one variable from active set. */
		case ST_INACTIVE:
			#ifndef __XPCTARGET__
			snprintf( messageString,80,"bound no. %d.", BC_idx );
			getGlobalMessageHandler( )->throwInfo( RET_REMOVE_FROM_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( removeBound( BC_idx,BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_REMOVE_FROM_ACTIVESET_FAILED );

			y[BC_idx] = 0.0;
			break;


		/* Add one variable to active set. */
		default:
			#ifndef __XPCTARGET__
			if ( BC_status == ST_LOWER )
				snprintf( messageString,80,"lower bound no. %d.", BC_idx );
			else
				snprintf( messageString,80,"upper bound no. %d.", BC_idx );
				getGlobalMessageHandler( )->throwInfo( RET_ADD_TO_ACTIVESET,messageString,__FUNCTION__,__FILE__,__LINE__,VS_VISIBLE );
			#endif

			if ( addBound( BC_idx,BC_status,BT_TRUE ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_ADD_TO_ACTIVESET_FAILED );
			break;
	}

	return SUCCESSFUL_RETURN;
}



/*
 * p e r f o r m D r i f t C o r r e c t i o n
 */
returnValue QProblemB::performDriftCorrection( )
{
	int i;
	int nV = getNV ();

	for ( i=0; i<nV; ++i )
	{
		switch ( bounds.getType ( i ) )
		{
			case ST_BOUNDED:
				switch ( bounds.getStatus ( i ) )
				{
					case ST_LOWER:
						lb[i] = x[i];
						ub[i] = getMax (ub[i], x[i]);
						y[i] = getMax (y[i], 0.0);
						break;
					case ST_UPPER:
						lb[i] = getMin (lb[i], x[i]);
						ub[i] = x[i];
						y[i] = getMin (y[i], 0.0);
						break;
					case ST_INACTIVE:
						lb[i] = getMin (lb[i], x[i]);
						ub[i] = getMax (ub[i], x[i]);
						y[i] = 0.0;
						break;
					case ST_UNDEFINED:
					case ST_INFEASIBLE_LOWER:
					case ST_INFEASIBLE_UPPER:
						break;
				}
				break;
			case ST_EQUALITY:
				lb[i] = x[i];
				ub[i] = x[i];
				break;
			case ST_UNBOUNDED:
			case ST_UNKNOWN:
            case ST_DISABLED:
				break;
		}
	}

	setupAuxiliaryQPgradient( );

	return SUCCESSFUL_RETURN;
}


/*
 *	s h a l l R e f a c t o r i s e
 */
BooleanType QProblemB::shallRefactorise( const Bounds* const guessedBounds ) const
{
	int i;
	int nV = getNV( );

	/* always refactorise if Hessian is not known to be positive definite */
	if ( getHessianType( ) == HST_SEMIDEF )
		return BT_TRUE;

	/* 1) Determine number of bounds that have same status
	 *    in guessed AND current bounds.*/
	int differenceNumber = 0;

	for( i=0; i<nV; ++i )
		if ( guessedBounds->getStatus( i ) != bounds.getStatus( i ) )
			++differenceNumber;

	/* 2) Decide wheter to refactorise or not. */
	if ( 2*differenceNumber > guessedBounds->getNFX( ) )
		return BT_TRUE;
	else
		return BT_FALSE;
}


/*
 *	a d d B o u n d
 */
returnValue QProblemB::addBound(	int number, SubjectToStatus B_status,
									BooleanType updateCholesky
									)
{
	int i, j;
	int nV  = getNV( );
	int nFR = getNFR( );


	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* Perform cholesky updates only if QProblemB has been initialised! */
	if ( getStatus( ) == QPS_PREPARINGAUXILIARYQP )
	{
		/* UPDATE INDICES */
		if ( bounds.moveFreeToFixed( number,B_status ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_ADDBOUND_FAILED );

		return SUCCESSFUL_RETURN;
	}


	/* I) PERFORM CHOLESKY UPDATE: */
	if ( ( updateCholesky == BT_TRUE ) &&
		 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
	{
		/* 1) Index of variable to be added within the list of free variables. */
		int number_idx = bounds.getFree( )->getIndex( number );

		real_t c, s, nu;

		/* 2) Use row-wise Givens rotations to restore upper triangular form of R. */
		for( i=number_idx+1; i<nFR; ++i )
		{
			computeGivens( RR(i-1,i),RR(i,i), RR(i-1,i),RR(i,i),c,s );
			nu = s/(1.0+c);

			for( j=(1+i); j<nFR; ++j ) /* last column of R is thrown away */
				applyGivens( c,s,nu,RR(i-1,j),RR(i,j), RR(i-1,j),RR(i,j) );
		}

		/* 3) Delete <number_idx>th column and ... */
		for( i=0; i<nFR-1; ++i )
			for( j=number_idx+1; j<nFR; ++j )
				RR(i,j-1) = RR(i,j);
		/* ... last column of R. */
		for( i=0; i<nFR; ++i )
			RR(i,nFR-1) = 0.0;
	}

	/* II) UPDATE INDICES */
	if ( bounds.moveFreeToFixed( number,B_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_ADDBOUND_FAILED );


	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e B o u n d
 */
returnValue QProblemB::removeBound(	int number,
									BooleanType updateCholesky
									)
{
	int i;
	int nV  = getNV( );
	int nFR = getNFR( );


	/* consistency check */
	if ( ( getStatus( ) == QPS_NOTINITIALISED )    ||
		 ( getStatus( ) == QPS_AUXILIARYQPSOLVED ) ||
		 ( getStatus( ) == QPS_HOMOTOPYQPSOLVED )  ||
		 ( getStatus( ) == QPS_SOLVED )            )
	{
		return THROWERROR( RET_UNKNOWN_BUG );
	}

	/* save index sets and decompositions for flipping bounds strategy */
	if ( options.enableFlippingBounds == BT_TRUE )
		flipper.set( &bounds,R );

	/* I) UPDATE INDICES */
	if ( bounds.moveFixedToFree( number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_REMOVEBOUND_FAILED );

	/* Perform cholesky updates only if QProblemB has been initialised! */
	if ( getStatus( ) == QPS_PREPARINGAUXILIARYQP )
		return SUCCESSFUL_RETURN;


	/* II) PERFORM CHOLESKY UPDATE */
	if ( ( updateCholesky == BT_TRUE ) &&
		 ( hessianType != HST_ZERO )   && ( hessianType != HST_IDENTITY ) )
	{
		int* FR_idx;
		bounds.getFree( )->getNumberArray( &FR_idx );

		/* 1) Calculate new column of cholesky decomposition. */
		real_t* rhs = new real_t[nFR+1];
		real_t* r   = new real_t[nFR];

		real_t r0;
		switch ( hessianType )
		{
			case HST_ZERO:
				if ( usingRegularisation( ) == BT_FALSE )
					r0 = 0.0;
				else
					r0 = options.epsRegularisation;
				for( i=0; i<nFR; ++i )
					rhs[i] = 0.0;
				break;

			case HST_IDENTITY:
				r0 = 1.0;
				for( i=0; i<nFR; ++i )
					rhs[i] = 0.0;
				break;

			default:
				H->getRow(number, bounds.getFree(), 1.0, rhs);
				r0 = H->diag(number);
				break;
		}

		if ( backsolveR( rhs,BT_TRUE,BT_TRUE,r ) != SUCCESSFUL_RETURN )
		{
			delete[] rhs; delete[] r;
			return THROWERROR( RET_REMOVEBOUND_FAILED );
		}

		for( i=0; i<nFR; ++i )
			r0 -= r[i]*r[i];

		/* 2) Store new column into R. */
		for( i=0; i<nFR; ++i )
			RR(i,nFR) = r[i];

		if ( options.enableFlippingBounds == BT_TRUE )
		{
			if ( r0 > options.epsFlipping )
				RR(nFR,nFR) = sqrt( r0 );
			else
			{
				hessianType = HST_SEMIDEF;

				flipper.get( &bounds,R );
				bounds.flipFixed(number);

				switch (bounds.getStatus(number))
				{
					case ST_LOWER: lb[number] = ub[number]; break;
					case ST_UPPER: ub[number] = lb[number]; break;
					default: delete[] rhs; delete[] r; return THROWERROR( RET_MOVING_BOUND_FAILED );
				}

			}
		}
		else
		{
			if ( r0 > ZERO )
				RR(nFR,nFR) = sqrt( r0 );
			else
			{
				delete[] rhs; delete[] r;

				hessianType = HST_SEMIDEF;
				return THROWERROR( RET_HESSIAN_NOT_SPD );
			}
		}

		delete[] rhs; delete[] r;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t I t e r a t i o n
 */
returnValue QProblemB::printIteration( 	int iteration,
										int BC_idx,	SubjectToStatus BC_status
		  								)
{
	#ifndef __XPCTARGET__
	char myPrintfString[160];

	/* consistency check */
	if ( iteration < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* nothing to do */
	if ( options.printLevel != PL_MEDIUM )
		return SUCCESSFUL_RETURN;


	/* 1) Print header at first iteration. */
 	if ( iteration == 0 )
	{
		snprintf( myPrintfString,160,"\n\n#################   qpOASES  --  QP NO. %3.0d   ##################\n\n", count );
		myPrintf( myPrintfString );

		myPrintf( "    Iter   |    StepLength    |       Info       |   nFX    \n" );
		myPrintf( " ----------+------------------+------------------+--------- \n" );
	}

	/* 2) Print iteration line. */
	if ( BC_status == ST_UNDEFINED )
	{
		if ( hessianType == HST_ZERO )
			snprintf( myPrintfString,80,"   %5.1d   |   %1.6e   |    LP SOLVED     |  %4.1d   \n", iteration,tau,getNFX( ) );
		else
			snprintf( myPrintfString,80,"   %5.1d   |   %1.6e   |    QP SOLVED     |  %4.1d   \n", iteration,tau,getNFX( ) );
		myPrintf( myPrintfString );
	}
	else
	{
		char info[8];

		if ( BC_status == ST_INACTIVE )
			snprintf( info,8,"REM BND" );
		else
			snprintf( info,8,"ADD BND" );

		snprintf( myPrintfString,80,"   %5.1d   |   %1.6e   |   %s %4.1d   |  %4.1d   \n", iteration,tau,info,BC_idx,getNFX( ) );
		myPrintf( myPrintfString );
	}
	#endif

	return SUCCESSFUL_RETURN;
}



END_NAMESPACE_QPOASES


/*
 *	end of file
 */
