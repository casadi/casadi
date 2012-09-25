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
 *	\file include/qpOASES/QProblem.ipp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0beta
 *	\date 2007-2012
 *
 *	Implementation of inlined member functions of the QProblem class which
 *	is able to use the newly developed online active set strategy for
 *	parametric quadratic programming.
 */


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	g e t C o n s t r a i n t s
 */
inline returnValue QProblem::getConstraints( Constraints& _constraints ) const
{
	int nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	_constraints = constraints;

	return SUCCESSFUL_RETURN;
}



/*
 *	g e t N C
 */
inline int QProblem::getNC( ) const
{
	return constraints.getNC( );
}


/*
 *	g e t N E C
 */
inline int QProblem::getNEC( ) const
{
	return constraints.getNEC( );
}


/*
 *	g e t N A C
 */
inline int QProblem::getNAC( ) const
{
	return constraints.getNAC( );
}


/*
 *	g e t N I A C
 */
inline int QProblem::getNIAC( ) const
{
	return constraints.getNIAC( );
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/


/*
 *	s e t A
 */
inline returnValue QProblem::setA( Matrix *A_new )
{
	int j;
	int nV = getNV( );
	int nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( A_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* Set constraint matrix AND update member AX. */
	if ( ( freeConstraintMatrix == BT_TRUE ) && ( A != 0 ) )
	{
		delete A;
		A = 0;
	}
	A = A_new;
	freeConstraintMatrix = BT_FALSE;

	A->times(1, 1.0, x, nV, 0.0, Ax, nC);

	for( j=0; j<nC; ++j )
	{
		Ax_u[j] = ubA[j] - Ax[j];
		Ax_l[j] = Ax[j] - lbA[j];
		
		// (ckirches) disable constraints with empty rows	
		if (A->getRowNorm (j) == 0.0)
			constraints.setType ( j, ST_DISABLED );
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	s e t A
 */
inline returnValue QProblem::setA( const real_t* const A_new )
{
	int j;
	int nV = getNV( );
	int nC = getNC( );
	DenseMatrix* dA;

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( A_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* Set constraint matrix AND update member AX. */
	if ( ( freeConstraintMatrix == BT_TRUE ) && ( A != 0 ) )
	{
		delete A;
		A = 0;
	}
	A = dA = new DenseMatrix(nC, nV, nV, (real_t*) A_new);
	freeConstraintMatrix = BT_TRUE;

	A->times(1, 1.0, x, nV, 0.0, Ax, nC);

	for( j=0; j<nC; ++j )
	{
		Ax_u[j] = ubA[j] - Ax[j];
		Ax_l[j] = Ax[j] - lbA[j];
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t L B A
 */
inline returnValue QProblem::setLBA( const real_t* const lbA_new )
{
	int i;
	int nV = getNV( );
	int nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( lbA_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	for( i=0; i<nC; ++i )
		lbA[i] = lbA_new[i];

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t L B A
 */
inline returnValue QProblem::setLBA( int number, real_t value )
{
	int nV = getNV( );
	int nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( ( number >= 0 ) && ( number < nC ) )
	{
		lbA[number] = value;
		return SUCCESSFUL_RETURN;
	}
	else
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );
}


/*
 *	s e t U B A
 */
inline returnValue QProblem::setUBA( const real_t* const ubA_new )
{
	int i;
	int nV = getNV( );
	int nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( ubA_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	for( i=0; i<nC; ++i )
		ubA[i] = ubA_new[i];

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t U B A
 */
inline returnValue QProblem::setUBA( int number, real_t value )
{
	int nV = getNV( );
	int nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( ( number >= 0 ) && ( number < nC ) )
	{
		ubA[number] = value;
		return SUCCESSFUL_RETURN;
	}
	else
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
