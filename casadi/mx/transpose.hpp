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

#ifndef TRANSPOSE_HPP
#define TRANSPOSE_HPP

#include "mx_node.hpp"
#include "reordering.hpp"

namespace CasADi{

/** Represents a transposition of an MX
  \author Joel Andersson 
  \date 2010
*/
class Transpose : public Reordering {
friend class MX;

public:

/** \brief  Constructor */
Transpose(const MX& x);

/** \brief  Clone function */
virtual Transpose* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

};

} // namespace CasADi

#endif // TRANSPOSE_HPP
