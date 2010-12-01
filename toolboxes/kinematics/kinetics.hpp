/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef KINETICS_HPP
#define KINETICS_HPP

#include "frame.hpp"
#include "kinvec.hpp"

namespace KINEMATICS{
using namespace CasADi;

/** \brief creates a 4x4 transformation matrix for a translation
* 
*   \param x The amount to shift over the x-axis
*   \param y The amount to shift over the y-axis
*   \param z The amount to shift over the z-axis
*
* \code
*  1 | 0 | 0 | x 
*  0 | 0 | 0 | y
*  0 | 0 | 1 | z
*  0 | 0 | 0 | 1
* \endcode
*/
SXMatrix translate(const SXMatrix &x,const SXMatrix &y,const SXMatrix &z);
/** \brief shorter notation for translate()
*/
SXMatrix tr(const SXMatrix &x,const SXMatrix &y,const SXMatrix &z);

/** \brief SXMatrix rotate(AXIS X,const SXMatrix &angle); // Rotate about an axis */
/** \brief creates a 3x3 rotation matrix for a rotation about the x-axis
* 
*   creates 3x3 rotation matrix for a rotation about the x-axis
*   \param angle angle for which to rotate in radians
*                Looking from the origin towards the endpoint of the x-defining unit vector, a positive rotation is clockwise
*
*   ca = cos(angle)
*   sa = sin(angle)
*
* \code
*  1 | 0  | 0  
*  0 | ca | -sa
*  0 | sa | ca  
* \endcode
*/
SXMatrix Rx(const SXMatrix &angle);
/** \brief creates a 3x3 rotation matrix for a rotation about the y-axis
* 
*   creates 3x3 rotation matrix for a rotation about the y-axis
*   \param angle angle for which to rotate in radians
*                Looking from the origin towards the endpoint of the y-defining unit vector, a positive rotation is clockwise
*
*   ca = cos(angle)
*   sa = sin(angle)
*
* \code
*  ca | -sa | 0 | 0 
*  sa | ca  | 0 | 0
*  0  | 0   | 0 | 0 
* \endcode
*/
SXMatrix Ry(const SXMatrix &angle);
/** \brief creates a 3x3 rotation matrix for a rotation about the z-axis
* 
*   creates 3x3 rotation matrix for a rotation about the z-axis
*   \param angle angle for which to rotate in radians
*                Looking from the origin towards the endpoint of the z-defining unit vector, a positive rotation is clockwise
*
*   ca = cos(angle)
*   sa = sin(angle)
*
* \code
*  1 | 0  | 0  
*  0 | ca | -sa
*  0 | sa | ca  
* \endcode
*/
SXMatrix Rz(const SXMatrix &angle);
/**
* \brief creates a 3x3 rotation matrix for a rotation about the x-axis with an angle that is a multiple of PI/2
* \param quadrant rotate over PI/2*quadrant
*/
SXMatrix Rxp(const int quadrant); // Rotate multiple of PI/2
/**
* \brief creates a 3x3 rotation matrix for a rotation about the y-axis with an angle that is a multiple of PI/2
* \param quadrant rotate over PI/2*quadrant
*/
SXMatrix Ryp(const int quadrant);
/**
* \brief creates a 3x3 rotation matrix for a rotation about the z-axis with an angle that is a multiple of PI/2
* \param quadrant rotate over PI/2*quadrant
*/
SXMatrix Rzp(const int quadrant);
SXMatrix TRx(const SXMatrix &angle);
SXMatrix TRy(const SXMatrix &angle);
SXMatrix TRz(const SXMatrix &angle);
SXMatrix TRxp(const int quadrant); // Rotate multiple of PI/2
SXMatrix TRyp(const int quadrant);
SXMatrix TRzp(const int quadrant);

/**
*  \brief Make a 3x3 rotation matrix that expresses a permutation of the axes
*  
* Make a 3x3 rotation matrix that expresses a permutation if the axes.
*
* Axes x,y,z are labeled as integers 1,2,3. A minus sign indicates a reversed direction
*
* \param a The new 1-axis is the old a-axis
* \param a The new 2-axis is the old b-axis
* \param a The new 3-axis is the old c-axis
*
*  The following example expresses a mirror operation around the z=0 plane:
* \code
*  Expression R=R(1,2,-3);
* \endcode
* Note that this would shift handedness of the frame. Not a good idea in mechanics. Better stick to conventional right-handed frames.
*/
SXMatrix Rperm(int a,int b,int c); // Rperm(a,b,c) permute axis. Axes x,y,z are labeled as 1,2,3
/** \brief  The new 1-axis is the old a-axis */
/** \brief  The new 2-axis is the old b-axis */
/** \brief  The new 3-axis is the old c-axis */
				     
				     /**
*  \brief Make a 4x4 transformation matrix that expresses a permutation of the axes
*
*  \see Rperm(int a,int b,int c)
*/
SXMatrix TRperm(int a,int b,int c);

} // namespace KINEMATICS


#endif // KINETICS_HPP
