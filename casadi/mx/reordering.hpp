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
  /** \brief Base class for different reorderings of matrices, such as transpose, slicing etc.

  We have:
  (l) input index
  (k) nonzero index of input
  (k') nonzero index of output 

  We map the nonzeros (k') of the output of the function to the input index (l) and the nonzero (k) of the input

  NOTE:
  Joel: Adapted to handle sparsity correctly, changed to vector for index mappings (for efficiency)

  \author Joris Gillis
  \date 2010
*/
class Reordering : public MXNode{
  public:

    /// Single input
    explicit Reordering(const MX &dep);

    /// Multiple inputs
    explicit Reordering(const std::vector<MX> &dep);

    /// Evaluate the function and store the result in the node
    virtual void evaluate(int fsens_order, int asens_order);

    /// Initialize
    virtual void init();

    /// Maps the non-zero of the result to the input index: Default all zeros
    int k2l(int k) const;

    /// Maps the non-zero of the result to the non-zero of the input
    int k2k(int k) const;

    /// Print
    virtual void print(std::ostream &stream=std::cout) const;

  protected:
    /// Mapping from the output non-zero to the input nonzero
    std::vector<int> nzind_;

    /// Mapping from the output non-zero to the input index
    std::vector<int> argind_;

};

} // namespace CasADi

#endif // REORDERING_HPP
