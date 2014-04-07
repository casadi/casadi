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

#ifndef PARALLELIZER_HPP
#define PARALLELIZER_HPP

#include <vector>

#include "function.hpp"

namespace CasADi{

  // Forward declaration of internal class
  class ParallelizerInternal;

  /** \brief Parallelizer execution of functions
      \author Joel Andersson
      \date 2011
  */ 
  class Parallelizer : public Function{
  public:

    /// Default constructor
    Parallelizer();

    /// Create a Parallelizer
    explicit Parallelizer(const std::vector<Function>& funcs);

    /// Access functions of the node
    ParallelizerInternal* operator->();

    /// Const access functions of the node
    const ParallelizerInternal* operator->() const;
  
  };

} // namespace CasADi


#endif // PARALLELIZER_HPP

