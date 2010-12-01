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

#ifndef IF_ELSE_NODE_HPP
#define IF_ELSE_NODE_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Represents a branch in an MX tree
  \author Joel Andersson 
  \date 2010
*/
class IfElseNode : public MXNode{
public:

/** \brief  Constructors */
IfElseNode(const MX& cond, const MX& if_true, const MX& if_false);
//IfElseNode(const MX& cond, const FX& fcn_true, const FX& fcn_false, const std::vector<MX>& arg);

/** \brief  Clone function */
virtual IfElseNode* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();
  
protected:

  // Functions to evaluate if true/false
//  FX fcn_true, fcn_false;
  
  //! \brief tolerance
  static double tol;
};

} // namespace CasADi


#endif // IF_ELSE_NODE_HPP
