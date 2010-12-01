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

#ifndef FRAME_HPP
#define FRAME_HPP
#include <casadi/expression_tools.hpp>
#include <string>
#include <iostream>

namespace KINEMATICS{
using namespace CasADi;

class FrameNode;
class KinVec;

/**
 *  \brief define the concept of a mechanical frame
 *  
 * A frame is a concept from classical mechanics.
 * One can define a frame {1} rigorously as an ordered tuple of four vectors: 
 *  - three 0-vectors (ex,ey,ez) denoting the axis orientations
 *  - one 1-vector (p) denoting the position of the frame's origin.
 *
 * It is convenient to write the constituents of {1} in a 4x4 matrix form T_10,
 * in which the defining vectors of frame {1} have been written in {0} coordinates:
 * T = [ex ey ez p]
 *
 *  or 
 * T = [ R  p ;0 0 0 1]
 *
 * with R a 3x3 rotation matrix
 *      p a 3x1 position vector
 *
 * In this toolbox, a Frame has a reference Frame and a Transformation matrix as members.
 * The user constructs a chain (or rather a tree) of frames by calling the constructors of Frame.
 *
 *
 * The longer explanation:
 *
 * Frames can be intuitively defined by the formalism of active interpretation of transformation matrices.
 * Put yourself in frame {i}. Now imagine taking a copy of the vector-tuplet that defines frame {i} and
 * moving and rotating these vectors at will. You have performed an action on objects in vector space.
 * You are at a certain pose of the copied tuple which you now say defines frame {j}.
 * 
 * You write down the vectors of {j} with respect to {i} and expressed in {i} in a matrix T_ji
 * It is convenient to use elemental transformations to obtain T_ji.
 * For example, {j} may be defined by starting from a copy of {i}, 
 * moving it over a distance L and rotating it over its (moved) z-axis:
 *
 * \code
 *	Expression T_ji=tr(L,0,0)*TRz(alpha);
 *	Expression T=TRz(alpha)*tr(L,0,0); // this is not the same transformation
 * \endcode
 *
 *  \author Joris Gillis
 *
 */
class Frame{
  friend class FrameNode;
  friend class KinVec;
  friend class LKinVec;
  public:
/** \brief  Constructors */
    Frame();
    /**
    * \brief Constructor for the world frame (top of the hierarchy)
    *
    * Constructor for the world frame (top of the hierarchy)
    * \param name	name of the Frame
    * \param q		nx1 expression containg the n default time-dependant symbols
    * \param dq		nx1 expression containg the derivatives of the n default time-dependant symbols
    * \param ddq	nx1 expression containg the second derivatives of the n default time-dependant variables
    * \param time	the symbol that represents time
    *
    * In the above, n is any integer, but the same at all three instances.
    */
    explicit Frame(const std::string& name,const SXMatrix &q,const SXMatrix &dq,const SXMatrix & ddq);
    
    /**
    * \brief Copy constructor
    * Copy constructor
    */
    Frame(Frame const& frame);

    /**
    * \brief Chain constructor for the world frame - to create the hierarchy of Frames
    *
    * Chain constructor for the world frame - to create the hierarchy of Frames
    * \param name	name of the Frame
    * \param ref	the reference Frame with respect to which the \a T holds
    * \param T		4x4 transformation matrix
    *			Can be conveniently constructed with the global functions TRx, TRy, ...
    **/
    Frame(const std::string& name, const Frame& ref,const SXMatrix & T); // T should be 4x4 homogenous matrix
    ~Frame();
    
    Frame& operator=(const Frame &frame);
    
    friend std::ostream& operator<<(std::ostream &stream, const Frame &frame);
    
    /** 
      \brief Get the name of this frame
    */
    const std::string& getName() const;
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
    Frame getCommonFrame(Frame other) const;  
    
    /** 
    *  \brief keep left-multiplying \a e by transformation matrices going from this frame to the \a ei frame 
    *
    * keep left-multiplying \a e by transformation matrices going from this frame to the \a ei frame 
    *  \param e		the expression that will be transformed
    *                   will typically be vector-like but can be a matrix as well
    *  \param ei	the Frame in which the result should be expressed (Expressed In)
    *  \param type	specifies if e is to be treated position-like (1) or direction-like (0)
    *  \return		the transformed expression
    *  \see KinVec
    *
    * Consider the Frame tree:  
    *\code
    *        ___0___
    *       1       2
    *     3   4       5
    * \endcode
    * (Frame 3).chain(e,(Frame 5))
    *
    * will return  inv(T_52)*inv(T_20)*T_10*T_31*e
    */
    SXMatrix chain(const SXMatrix &e,const Frame &ei,bool type) const;
  protected:
  FrameNode* node;
    explicit Frame(FrameNode* ptr);
    
    std::set<FrameNode*> getReferences() const;	// Get chain of references
    FrameNode* getCommonFramePtr(Frame other) const;  // Get the nearest ancestor that connects 2 frames
    std::vector<FrameNode*> getFrameChain(FrameNode* endptr) const;  // Bubble upwards until endpoint is met
    

};

} // namespace KINEMATICS


#endif //FRAME_HPP

