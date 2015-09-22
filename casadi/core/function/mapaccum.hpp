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


#ifndef CASADI_MAPACCUM_HPP
#define CASADI_MAPACCUM_HPP

#include "function.hpp"

namespace casadi {

  /** \brief  Forward declaration of internal class */
  class MapAccumInternal;

  /** \brief MapAccum class

      Consider a function: f(x, u) -> xp, y

      MapAccum will evaluate this function repeatedly,
      feeding the output back to the input:

      \verbatim
        x1, y0 <- f(x0, u0)
        x2, y1 <- f(x1, u1)
        x3, y2 <- f(x2, u2)
      \endverbatim

      The inputs to MapAccum(f) are in this case:
        - x0
        - [u0 u1 u2]
      The outputs are:
        - [x1 x2 x3]
        - [y0 y1 y2]

      This class treats the general case:
        - n repetitions
        - any number of accumulated inputs (x)
        - any number of regular inputs (u)
        - a boolean list `input_accum` flags which inputs are accumalated.
        - An index list `output_accum` indicates the indices of the outputs
          that are fed back to the inputs.

      This implementation is optimized for speed.
      There is a penalty in memory:
        the entire accumulator history is an output.
      This allows the forward mode to use this history,
      instead of recreating the accumulator.

      In reverse mode, you would need the history anyway.

      For a memory-optimized implementation, see Fold.

      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT MapAccum : public Function {
  public:
    /** \brief Default constructor */
    MapAccum();

    /** \brief Constructor (generic mapaccum) */
    MapAccum(const std::string& name, const Function& f,
           int n,
           const std::vector<bool>& input_accum,
           const std::vector<int>& output_accum,
           bool reverse = false,
           const Dict& opts=Dict());


    /** \brief  Access functions of the node */
    MapAccumInternal* operator->();

    /** \brief  Const access functions of the node */
    const MapAccumInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };

} // namespace casadi

#endif // CASADI_MAPACCUM_HPP
