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

#ifndef C_FUNCTION_HPP
#define C_FUNCTION_HPP

#include "fx.hpp"
#include <string>

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class CFunctionNode;

// Forward declaration
class CFunction;

/** \brief  Wrapper around functions */
typedef void (*CFunctionWrapper)(CFunction &f, int fsens_order, int asens_order, void* user_data);

/** \brief  Interface to function implemented as plan code 
  \author Joel Andersson 
  \date 2010
*/
class CFunction : public FX{

public:

/** \brief  default constructor */
  CFunction();

  /** \brief  Create a function */
  explicit CFunction(CFunctionWrapper c_fcn);

  /** \brief  Access functions of the node */
  CFunctionNode* operator->();
  const CFunctionNode* operator->() const;
  
  /** \brief  Make sure that the pointer points towards a valid object */
  void assertNode() const;
  
}; // class CFunction

/** \brief  Internal class for CFunction
  \author Joel Andersson 
  \date 2010
 A regular user should never work with any Node class. Use CFunction directly.
*/
class CFunctionNode : public FXNode{
  friend class CFunction;
  public:
    
    /** \brief  Create a function */
    explicit CFunctionNode(CFunctionWrapper c_fcn);
    
    /** \brief  Destructor */
    virtual ~CFunctionNode();

    /** \brief  Set user data structure (to be passed to all functions) */
    void setUserData(void* user_data);
  
    /** \brief  Evaluate */
    virtual void evaluate(int fsens_order, int asens_order);
  
    /** \brief  Initialize */
    virtual void init();
  
  
  protected:
    void* user_data_;
    CFunctionWrapper evaluate_;
  
    /// A reference to this object to be passed to the user functions
    CFunction ref_;
    
}; // class CFunctionNode 
  


} // namespace CasADi


#endif // C_FUNCTION_HPP
