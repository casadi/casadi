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
 *	\file include/qpOASES/Utils.ipp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0beta
 *	\date 2007-2012
 *
 *	Implementation of some inlined utilities for working with the different QProblem classes.
 */


BEGIN_NAMESPACE_QPOASES


/*
 *   g e t S i g n
 */
inline real_t getSign(	real_t arg
						)
{
	if ( arg >= 0.0 )
		return 1.0;
	else
		return -1.0;
}



/*
 *   g e t M a x
 */
inline int getMax(	int x,
					int y
					)
{
    return (y<x) ? x : y;
}


/*
 *   g e t M i n
 */
inline int getMin(	int x,
					int y
					)
{
    return (y>x) ? x : y;
}



/*
 *   g e t M a x
 */
inline real_t getMax(	real_t x,
						real_t y
						)
{
    return (y<x) ? x : y;
}


/*
 *   g e t M i n
 */
inline real_t getMin(	real_t x,
						real_t y
						)
{
    return (y>x) ? x : y;
}


/*
 *   g e t A b s
 */
inline real_t getAbs(	real_t x
						)
{
    return (x>=0) ? x : -x;
}



END_NAMESPACE_QPOASES


/*
 *	end of file
 */
