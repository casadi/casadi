/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef SYMBOLIC_SXElement_HPP
#define SYMBOLIC_SXElement_HPP

#include "sx_node.hpp"
/// \cond INTERNAL

namespace casadi {

/** \brief Represents a scalar symbolic expression
  \author Joel Andersson
  \date 2010
  A regular user is not supposed to work with this Node class.
  This user can call SX(name) instead.
*/
class CASADI_CORE_EXPORT SymbolicSX : public SXNode {
public:
  explicit SymbolicSX(const std::string &name) : name(name) {}
  virtual ~SymbolicSX() {}

  virtual bool isSymbolic() const { return true; }

  virtual const std::string& getName() const { return name; }

    /** \brief  Get the operation */
  virtual int getOp() const { return OP_PARAMETER;}

  /** \brief  Name */
  std::string name;

protected:

  /** \brief  print */
  virtual void print(std::ostream &stream, long& remaining_calls) const {
    stream << name;
  }
};

} // namespace casadi
/// \endcond
#endif // SYMBOLIC_SXElement_HPP
