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

#ifndef SLICE_HPP
#define SLICE_HPP

#include "mx_node.hpp"
#include "reordering.hpp"
#include "slicer.hpp"

namespace CasADi{

/** Represents a slice of an MX
  \author Joris Gillis 
  \date 2010
*/
class Slice : public Reordering {
friend class MX;

public:

/** \brief  Constructor */
Slice(const MX& x, Slicer i, Slicer j);

/** \brief  Clone function */
virtual Slice* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/// Initialize
virtual void init();

/// Evaluate (do you really need this?)
virtual void evaluate(int fsens_order, int asens_order);

protected:
Slicer i;
Slicer j;
};

} // namespace CasADi

#endif // SLICE_HPP
