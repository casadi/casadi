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

#ifndef REORDERING_HPP
#define REORDERING_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Base class for MXNodes that are a mere reoderening of the entries of their dependencies.

  This class mainly deals with indexing magic.
  
  We have:
  (i,j) matrix type indexing.
  (k) vector type indexing.
  
  We map k to a dependency (l*) and a vector index for this dependency (l)

  \author Joris Gillis
  \date 2010
*/
class Reordering : public MXNode{
public:

/** \brief  Constructor */
explicit Reordering(const std::vector<MX> &comp);

explicit Reordering(const MX &comp);

/** \brief  Evaluate the function and store the result in the node
 There is a default implementation that works with the indexing functions.
 It will work, but for speed, reimplement method.
 */
virtual void evaluate(int fsens_order, int asens_order);

/** \brief Maps (k)  to (l)
Default implementation: return 0
*/
virtual int k2l(int k) { return 0;};

/** \brief Maps (k)  to (k*)
*/
virtual int k2k(int k)=0;

// These index remapping function out to be inlined. Can one inline without losing polymorphism here?
// Some C++ guru ought to take a look here...

protected:
  
  
};

} // namespace CasADi

#endif // REORDERING_HPP
