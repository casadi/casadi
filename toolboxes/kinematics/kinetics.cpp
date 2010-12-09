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

#include "kinetics.hpp"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>


enum AXIS{AX_X,AX_Y,AX_Z};

namespace KINEMATICS{
using namespace std;
using namespace CasADi;

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
SXMatrix Rperm(int a, int b, int c) {
  SXMatrix R(3,3);
  R(abs(a)-1,0)=a/abs(a);
  R(abs(b)-1,1)=b/abs(b);
  R(abs(c)-1,2)=c/abs(c);
  return R;
}

SX cosquadrant(const int quadrant) {
  int k=quadrant % 4; 
  if (k < 0) k+=4;// force mod to be positive
  SX e=SX::zero;
  if (k==0) e=SX::one;
  if (k==2) e=-SX::one;
  return e;
}

SX sinquadrant(const int quadrant) {
  int k=quadrant % 4; 
  if (k < 0) k+=4;// force mod to be positive
  SX e=SX::zero;
  if (k==1) e=SX::one;
  if (k==3) e=-SX::one;
  return e;
}

/**
* \brief creates a 3x3 rotation matrix for a rotation about the x-axis with an angle that is a multiple of PI/2
* \param quadrant rotate over PI/2*quadrant
*/
SXMatrix Rxp(const int quadrant) {
  SXMatrix R(3,3);
  SXMatrix ca=cosquadrant(quadrant);SXMatrix sa=sinquadrant(quadrant);
  //matrix([1,0,0],[0,cos(alpha),-sin(alpha)],[0,sin(alpha),cos(alpha)])
  R(0,0)=SX::one;R(1,1)=ca;R(1,2)=-sa;R(2,1)=sa;R(2,2)=ca;
  return R;
}

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
SXMatrix Rx(const SXMatrix &angle) {
  SXMatrix R(3,3);
  SXMatrix ca=cos(angle);SXMatrix sa=sin(angle);
  //matrix([1,0,0],[0,cos(alpha),-sin(alpha)],[0,sin(alpha),cos(alpha)])
  R(0,0)=SX::one;R(1,1)=ca;R(1,2)=-sa;R(2,1)=sa;R(2,2)=ca;
  return R;
}

/**
* \brief creates a 3x3 rotation matrix for a rotation about the y-axis with an angle that is a multiple of PI/2
* \param quadrant rotate over PI/2*quadrant
*/
SXMatrix Ryp(const int quadrant) {
  SXMatrix R(3,3);
  SXMatrix ca=cosquadrant(quadrant);SXMatrix sa=sinquadrant(quadrant);
  //matrix([cos(alpha),0,sin(alpha)],[0,1,0],[-sin(alpha),0,cos(alpha)])
  R(0,0)=ca;R(0,2)=sa;R(1,1)=SX::one;R(2,0)=-sa;R(2,2)=ca;
  return R;
}
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
SXMatrix Ry(const SXMatrix &angle) {
  SXMatrix R(3,3);
  SXMatrix ca=cos(angle);SXMatrix sa=sin(angle);
  //matrix([cos(alpha),0,sin(alpha)],[0,1,0],[-sin(alpha),0,cos(alpha)])
  R(0,0)=ca;R(0,2)=sa;R(1,1)=SX::one;R(2,0)=-sa;R(2,2)=ca;
  return R;
}

/**
* \brief creates a 3x3 rotation matrix for a rotation about the z-axis with an angle that is a multiple of PI/2
* \param quadrant rotate over PI/2*quadrant
*/
SXMatrix Rzp(const int quadrant) {
  SXMatrix R(3,3);
  SXMatrix ca=cosquadrant(quadrant);SXMatrix sa=sinquadrant(quadrant);
  //matrix([cos(alpha),0,sin(alpha)],[0,1,0],[-sin(alpha),0,cos(alpha)])
  R(0,0)=ca;R(0,1)=-sa;R(1,0)=sa;R(1,1)=ca;R(2,2)=SX::one;
  return R;
}
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
SXMatrix Rz(const SXMatrix &angle) {
  SXMatrix R(3,3);
  SXMatrix ca=cos(angle);SXMatrix sa=sin(angle);
  //matrix([cos(alpha),0,sin(alpha)],[0,1,0],[-sin(alpha),0,cos(alpha)])
  R(0,0)=ca;R(0,1)=-sa;R(1,0)=sa;R(1,1)=ca;R(2,2)=SX::one;
  return R;
}

SXMatrix TRx(const SXMatrix &angle) {
  SXMatrix T(4,4);
  T(3,3)=SX::one;
  SXMatrix R=Rx(angle);
  T(0,0)=R(0,0);T(0,1)=R(0,1);T(0,2)=R(0,2);
  T(1,0)=R(1,0);T(1,1)=R(1,1);T(1,2)=R(1,2);
  T(2,0)=R(2,0);T(2,1)=R(2,1);T(2,2)=R(2,2);
  return T;
}

SXMatrix TRy(const SXMatrix &angle) {
  SXMatrix T(4,4);
  T(3,3)=SX::one;
  SXMatrix R=Ry(angle);
  T(0,0)=R(0,0);T(0,1)=R(0,1);T(0,2)=R(0,2);
  T(1,0)=R(1,0);T(1,1)=R(1,1);T(1,2)=R(1,2);
  T(2,0)=R(2,0);T(2,1)=R(2,1);T(2,2)=R(2,2);
  return T;
}

SXMatrix TRz(const SXMatrix &angle) {
  SXMatrix T(4,4);
  T(3,3)=SX::one;
  SXMatrix R=Rz(angle);
  T(0,0)=R(0,0);T(0,1)=R(0,1);T(0,2)=R(0,2);
  T(1,0)=R(1,0);T(1,1)=R(1,1);T(1,2)=R(1,2);
  T(2,0)=R(2,0);T(2,1)=R(2,1);T(2,2)=R(2,2);
  return T;
}

SXMatrix TRxp(const int quadrant) {
  SXMatrix T(4,4);
  T(3,3)=SX::one;
  SXMatrix R=Rxp(quadrant);
  T(0,0)=R(0,0);T(0,1)=R(0,1);T(0,2)=R(0,2);
  T(1,0)=R(1,0);T(1,1)=R(1,1);T(1,2)=R(1,2);
  T(2,0)=R(2,0);T(2,1)=R(2,1);T(2,2)=R(2,2);
  return T;
}

SXMatrix TRyp(const int quadrant) {
  SXMatrix T(4,4);
  T(3,3)=SX::one;
  SXMatrix R=Ryp(quadrant);
  T(0,0)=R(0,0);T(0,1)=R(0,1);T(0,2)=R(0,2);
  T(1,0)=R(1,0);T(1,1)=R(1,1);T(1,2)=R(1,2);
  T(2,0)=R(2,0);T(2,1)=R(2,1);T(2,2)=R(2,2);
  return T;
}

SXMatrix TRzp(const int quadrant) {
  SXMatrix T(4,4);
  T(3,3)=SX::one;
  SXMatrix R=Rzp(quadrant);
  T(0,0)=R(0,0);T(0,1)=R(0,1);T(0,2)=R(0,2);
  T(1,0)=R(1,0);T(1,1)=R(1,1);T(1,2)=R(1,2);
  T(2,0)=R(2,0);T(2,1)=R(2,1);T(2,2)=R(2,2);
  return T;
}
/**
*  \brief Make a 4x4 transformation matrix that expresses a permutation of the axes
*
*  \see Rperm(int a,int b,int c)
*/
SXMatrix TRperm(int a, int b, int c) {
  SXMatrix T(4,4);
  T(3,3)=SX::one;
  SXMatrix R=Rperm(a,b,c);
  T(0,0)=R(0,0);T(0,1)=R(0,1);T(0,2)=R(0,2);
  T(1,0)=R(1,0);T(1,1)=R(1,1);T(1,2)=R(1,2);
  T(2,0)=R(2,0);T(2,1)=R(2,1);T(2,2)=R(2,2);
  return T;
}

/** \brief shorter notation for translate()
*/
SXMatrix tr(const SXMatrix &x,const SXMatrix &y,const SXMatrix &z) {
  return translate(x,y,z);
}

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
SXMatrix translate(const SXMatrix &x,const SXMatrix &y,const SXMatrix &z) {
  SXMatrix T(4,4);
  T(0,0)=SX::one;
  T(1,1)=SX::one;
  T(2,2)=SX::one;
  T(0,3)=x;T(1,3)=y;T(2,3)=z;
  T(3,3)=SX::one;
  return T;
}

// SXMatrix rotate (AXIS X,const SXMatrix &angle) {
//   switch ( X ) {
//     case AX_X : 
//       return TRx(angle);
//       break;
//     case AX_Y : 
//       return TRy(angle);
//       break;
//     case AX_Z : 
//       return TRz(angle);
//       break;
//   }
//  
// }


} // namespace KINEMATICS
