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

#ifndef SYMBOLIC_MATRIX_HPP
#define SYMBOLIC_MATRIX_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents a symbolic MX
  \author Joel Andersson 
  \date 2010
  A regular user is not supposed to work with this Node class.
  This user can call MX(name,n,m) directly.
*/
class SymbolicMatrix : public MXNode{
public:

/** \brief  Constructors */
explicit SymbolicMatrix(const std::string& name, int n=1, int m=1);

/** \brief  Clone function */
virtual SymbolicMatrix* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();

/** \brief  Is symbolic */
  virtual bool isSymbolic() const;
  
/** \brief  Get the name */
  virtual const std::string& getName() const;
  
protected:
  std::string name;
};

} // namespace CasADi


#endif // SYMBOLIC_MATRIX_HPP
