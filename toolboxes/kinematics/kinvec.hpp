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

#ifndef KINVEC_HPP
#define KINVEC_HPP
#include <casadi/expression_tools.hpp>
#include <string>
#include <iostream>
#include "frame.hpp"
#include "kinvec.hpp"

namespace KINEMATICS{
using namespace CasADi;

/**
 *  \brief Represent kinematical vectors
 *
 * What is a kinematic vector?
 *
 * A kinematic vector is a - possibly time varying - object in vectorspace that can exist irrespective of what
 * reference frames happen to be defined in this space. In other words, a vector has an identity of its own
 * and a given column matrix r=[r_x;r_y;r_z] (or vector as it is often confusingly called) is merely a image
 * of this vector from one point of view (expressed or resolved or represented in one particular frame). 
 *
 * In fact two types of vectors can be distinguished: position or 1-vectors and displacement or 0-vectors.
 * - For a 0-vector, only orientation and length matter.
 * - For a 1-vector, the anchor point of the vector is important part of the vector's definition.
 *
 * This class represents a kinematic vector by having a column matrix representation and a reference frame as members.
 *
 * Operators can work on KinVec objects with different reference frames. The result is automatically computed in a common frame.
 *
 * The really useful functions are found in the header documentation (KinVec.hpp)
 *
 *
 * \code
 *   KinVec a=rotVel(f1,f0,f0)+rotVel(f2,f1,f0);
 *   KinVec b=rotVel(f2,f0,f0); // a-c is zero (relation from mechanics)
 *   KinVec c=rotVel(f2,f0,f1); // b-c is zero
 * \endcode
 * 
 *  \author Joris Gillis
 *
 */
class KinVec {
  public:
    KinVec();
    KinVec(const KinVec &v);
    KinVec(const SXMatrix& x,const SXMatrix& y,const SXMatrix& z, bool type, const Frame& ref);
    KinVec(const SXMatrix& xyz,bool type,const Frame& ref);
    KinVec(const SXMatrix& xyz, bool type, const Frame& ref,const SXMatrix& J,const SXMatrix & q,const SXMatrix & dq,const SXMatrix & ddq,int order);
    KinVec(const SXMatrix& xyz, bool type, const Frame& ref,const SXMatrix& J,const SXMatrix &c,const SXMatrix & q,const SXMatrix & dq,const SXMatrix & ddq,int order);
    KinVec(const SXMatrix& xyz, bool type, const Frame& ref,const SXMatrix & q,const SXMatrix & dq,const SXMatrix & ddq,int order);
    
    
    friend std::ostream& operator<<(std::ostream &stream, const KinVec &vec);


    /** \brief Return the components of the KinVec as 3x1 expression
    *
    * \return 3x1 Expression
    */
    SXMatrix getCoords() const; // Get coordinates as an SXMatrix(3)
    /**
    * \return 0 for velocity vectors, one for position vectors
    */
    bool getType(); // 0 for velocity vectors, one for position vectors
    Frame getRef();
    
    /** 
    *  \brief Get the default time dependant symbols from the world frame
    */
    const SXMatrix& getQ() const;
    /** 
    *  \brief Get the derivatives of the default time dependant symbols from the world frame
    */
    const SXMatrix& getDQ() const;
    /** 
    *  \brief Get the second derivatives of the default time dependant symbols from the world frame
    */
    const SXMatrix& getDDQ() const;
    /**  
    *  \brief returns the lowermost Frame in the hierarchy this Frame shares with the \a other Frame
    */
    
    void setDDQ(const SXMatrix& ddq_);
    
    friend KinVec operator+(const KinVec &a,const KinVec &b);
    friend SXMatrix operator*(const KinVec &a,const KinVec &b);
    friend KinVec operator-(const KinVec &a,const KinVec &b);
    /** \brief Take the cross product of two vectors */
    friend KinVec cross(const KinVec &a,const KinVec &b);
    
    friend KinVec operator*(const KinVec &a,const SXMatrix &b);
    friend KinVec operator*(const SXMatrix &a,const KinVec &b);
    friend KinVec operator/(const KinVec &a,const SXMatrix &b);
    
    KinVec operator-() const;
    
    KinVec& operator+=(const KinVec &b);
    KinVec& operator-=(const KinVec &b);
    KinVec operator-(const KinVec &b);
	
    /** \brief Take the 2-norm of a vector */
    friend SXMatrix norm(const KinVec &v);
    void splitdep(const SXMatrix &x,const SXMatrix &dx,KinVec &v1,KinVec &v2 );
    /** \brief express the KinVec in another Frame */
    KinVec expressedIn(const Frame& f);
    /** \brief component-wise time derivative */
    KinVec der(); // component-wise derivative
        /**
    \param ddq symbols for the second derivatives, to which is explicitized
    */
    SXMatrix explicitize(const SXMatrix & q,const SXMatrix & ddq) const; // ddq contains symbols for the second derivatives, it should contain exaclty three nans, which signify the variables to which to explicitize
    SXMatrix explicitize(const SXMatrix & qi) const;
    
        /**
    \param ddq symbols for the second derivatives, to which is explicitized
    */
    SXMatrix explicitize(std::map<int,int> &di) const;
    SXMatrix explicitize(std::vector<int> &di) const;
    
    /**
    * \param di map which contains indices of ddq to which is explicitized
    * \param ri map which contains indices of the rows of the KinVec which are taken into account
    *
    * Documentation not finished
    * 
    * Will construct a system A*x+b from KinVec
    * eg:
    * vec:
    * \code
    * | 1 0 0 0 | | ddr     |   | x |
    * | 0 1 0 0 | | ddphi   | + | y |
    * | 1 1 1 1 | | ddtheta |   | z |
    *             | dddelta |
    * \endcode
    *
    * ddqn=[0 0 0 0];
    */
    SXMatrix explicitize(std::map<int,int> &di,std::map<int,int> &ri) const;
    SXMatrix explicitize(std::vector<int> &di,std::vector<int> &ri) const;
    
    
    /** \brief vector style access (allocation if element does not exist)
    */
    SXMatrix::Element operator[](int i);                      // vector style access (allocation if element does not exist)
    /** \brief vector style access (never allocate) 
    */
    const SX operator[](int i) const;          // vector style access (never allocate)

    SXMatrix v;
    SXMatrix q;
    SXMatrix dq;
    SXMatrix ddq;
    SXMatrix J;
    SXMatrix c;
    int order;
    bool type; // 0 for velocity vectors, one for position vectors
    Frame ref;
    
  protected:
    /** \brief Take two vectors, express them in a common frame and return this common frame.
    *
    * Take two vectors, express them in a common frame and return this common frame.
    *
    * Will change the KinVecs passed to it
    */
    friend Frame expressCommon(KinVec &a,KinVec &b);


};


/** \brief  Create a kinvec */
/** \brief  Get properties of the origin */
/**
* \brief  Get the 1-vector that defines the position a frame \a f expressed it in another frame \a ei
*
* Get the 1-vector that defines the position a frame \a f expressed it in another frame \a ei
*
* Example usage:
*
* \code
*   Frame f1("f1",f0,tr(x,y,z));
*   KinVec p=pos(f1,f1);        // Would be [0,0,0,1]
*   KinVec q=p.ExpressedIn(f0); // Would be [x,y,z,1]
*   KinVec r=pos(f1,f0);	// same as above
* \endcode
*/
KinVec pos(const Frame& f,const Frame& ei);		// position
/**
* \brief  Get the linear velocity of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* Get the velocity of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
* \param q		nx1 expression containg the n time-dependant symbols
* \param dq		nx1 expression containg the derivatives of the n time-dependant symbols
* \param ddq		nx1 expression containg the second derivatives of the n time-dependant variables
*
*/
KinVec vel(const Frame& f,const Frame& wt, const Frame& ei,const SXMatrix & q,const SXMatrix & dq,const SXMatrix & ddq);     // linear velocity
/**
* \brief  Get the linear velocity of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei using default time-dependence
*
* Get the velocity of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* The default time-dependant symbols are taking from the world frame
*
* Example usage:
*
* \code
*   Frame f1("f1",f0,TRz(a)*tr(0,0,H));
*   Frame f2("f2",f1,tr(L,0,0));
*   KinVec	p=vel(f2,f1,f1);      	// Would be [dL,0,0,0]
*   p=pos(f2,f1,f1).der(t); 		// Same
*
*   p=vel(f2,f0,f1);      		// Would be [dL,da*L,0,0]
*   p=vel(f2,f0,f0).ExpressedIn(f1); 	// Same
*
*   p=vel(f2,f1,f0);    		// Would be [cos(a)*dL,sin(a)*dL,0,0]
*   p=vel(f2,f1,f1).ExpressedIn(f0); 	// Same
*
*   p=vel(f2,f0,f0);     		// Would be [cos(a)*dL-sin(a)*da*L,sin(a)*dL+cos(a)*da*L,0,0]
*   p=pos(f2,f0,f0).der(t); 		// Same
* \endcode
*
* Note that - in general - velocity is not just the component-wise time derivative of positition.
* This only holds if \a ei equals \a wt
*
*/
KinVec vel(const Frame& f,const Frame& wt, const Frame& ei);    
/**
* \brief  Get the linear acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* Get the acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
* \param q		nx1 expression containg the n time-dependant symbols
* \param dq		nx1 expression containg the derivatives of the n time-dependant symbols
* \param ddq		nx1 expression containg the second derivatives of the n time-dependant variables
*/
KinVec acc(const Frame& f,const Frame& wt, const Frame& ei,const SXMatrix & q,const SXMatrix & dq);     // linear acceleration
/**
* \brief  Get the linear acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei using default time-dependence
*
* Get the linear acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* The default time-dependant symbols are taking from the world frame
*
* \code
*   KinVec	p=acc(f2,f1,f1);      	
*   p=vel(f2,f1,f1).der(t); 		// Same as previous line
*
*   p=acc(f2,f0,f1);      		
*   p=acc(f2,f0,f0).ExpressedIn(f1); 	// Same as previous line
*
*   p=acc(f2,f1,f0);    		 
*   p=acc(f2,f1,f1).ExpressedIn(f0); 	// Same as previous line
*
*   p=acc(f2,f0,f0);     		
*   p=vel(f2,f0,f0).der(t); 		// Same as previous line
* \endcode
*
* Note that - in general - acceleration is not just the component-wise time derivative of velocity.
* This only holds if \a ei equals \a wt
*/
KinVec acc(const Frame& f,const Frame& wt, const Frame& ei);
/**
* \brief  Get the rotational velocity of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* Get the rotational velocity of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
* \param q		nx1 expression containg the n time-dependant symbols
* \param dq		nx1 expression containg the derivatives of the n time-dependant symbols
* \param ddq		nx1 expression containg the second derivatives of the n time-dependant variables
*/
KinVec rotVel(const Frame& f,const Frame& wt,const  Frame& ei,const SXMatrix & q,const SXMatrix & dq);   // rotational velocity
/**
* \brief  Get the rotational velocity of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei using default time-dependence
*
* Get the rotational of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* The default time-dependant symbols are taking from the world frame.
*
* The rotational velocity is defined by its skew-symmetric matrix form:
*    W_10 = dot(R_10)*R_01
* 
*    with R_10 the 3x3 rotation matrix that transforms from frame 0 to 1.
*/

KinVec rotVel(const Frame& f,const Frame& wt,const  Frame& ei); 
/**
* \brief  Get the rotational acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* Get the rotational acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
* \param q		nx1 expression containg the n time-dependant symbols
* \param dq		nx1 expression containg the derivatives of the n time-dependant symbols
* \param ddq		nx1 expression containg the second derivatives of the n time-dependant variables
*/
KinVec rotAcc(const Frame& f,const Frame& wt,const  Frame& ei,const SXMatrix & q,const SXMatrix & dq);	// rotational acceleration
/**
* \brief  Get the rotational acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei using default time-dependence
*
* Get the rotational acceleration of the origin of frame \a f with respect to frame \a wt, expressed in frame \a ei
*
* The default time-dependant symbols are taking from the world frame
*
* Note that - contrary to the case with linear acceleration - rotational acceleration IS just the component-wise time derivative of rotational velocity.
*/
KinVec rotAcc(const Frame& f,const Frame& wt,const  Frame& ei);

/**
* \brief  Get the 0-vector that defines the x-axis of a frame \a f
*
* Get the 0-vector that defines the x-axis of a frame \a f
*/
KinVec ex(const Frame& f);
/**
* \brief  Get the 0-vector that defines the y-axis of a frame \a f
*
* Get the 0-vector that defines the x-axis of a frame \a f
*/
KinVec ey(const Frame& f);
/**
* \brief  Get the 0-vector that defines the z-axis of a frame \a f
*
* Get the 0-vector that defines the x-axis of a frame \a f
*/
KinVec ez(const Frame& f);
/**
* \brief  Get the 0-vector that defines the x-axis of a frame \a f, but express it in another frame \a ei
*
* Get the 0-vector that defines the x-axis of a frame \a f, but express it in another frame \a ei
* \code
* KinVec e0=ex(f1).expressedIn(ei);
* KinVec e1=ex(f1,ei); // same as e0
* \endcode
*/
KinVec ex(const Frame& f,const Frame& ei);
/**
* \brief  Get the 0-vector that defines the y-axis of a frame \a f, but express it in another frame \a ei
*
* Get the 0-vector that defines the y-axis of a frame \a f, but express it in another frame \a ei
* \code
* KinVec e0=ey(f1).expressedIn(ei);
* KinVec e1=ey(f1,ei); // same as e0
* \endcode
*/
KinVec ey(const Frame& f,const Frame& ei);
/**
* \brief  Get the 0-vector that defines the z-axis of a frame \a f, but express it in another frame \a ei
*
* Get the 0-vector that defines the z-axis of a frame \a f, but express it in another frame \a ei
* \code
* KinVec e0=ez(f1).expressedIn(ei);
* KinVec e1=ez(f1,ei); // same as e0
* \endcode
*/
KinVec ez(const Frame& f,const Frame& ei);


} // namespace KINEMATICS

#endif //KINVEC_HPP

