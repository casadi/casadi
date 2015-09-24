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


#ifndef CASADI_MAP_HPP
#define CASADI_MAP_HPP

#include "function.hpp"

namespace casadi {

  /** \brief  Forward declaration of internal class */
  class MapReduce;

  /** Map 
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT Map : public Function {
  public:
    /** \brief Default constructor */
    Map();

    /** \brief Constructor (any type of map) */
    Map(const std::string& name, const Function& f, int n,
        const Dict& opts=Dict());

    /** \brief Constructor (generic map) */
    Map(const std::string& name, const Function& f,
           int n,
           const std::vector<bool> &repeat_in,
           const std::vector<bool> &repeat_out,
           const Dict& opts=Dict());

    /** \brief  Access functions of the node */
    MapReduce* operator->();

    /** \brief  Const access functions of the node */
    const MapReduce* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };

} // namespace casadi

#endif // CASADI_MAP_HPP
