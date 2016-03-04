
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


#ifndef CASADI_JIT_HPP
#define CASADI_JIT_HPP

#include "function.hpp"
#include "../casadi_file.hpp"

namespace casadi {

  /** \brief Create a just-in-time compiled function from a C/C++ language string
   * The function can an arbitrary number of inputs and outputs that must all be
   * scalar-valued.
   * Only specify the function body, assuming that the inputs are stored in an array
   * named 'arg' and the outputs stored in an array named 'res'. The data type
   * used must be 'real_t', which is typically equal to 'double` or another data
   * type with the same API as 'double'.
   *
   * The final generated function will have a structure similar to:
   * 
   * void fname(const real_t* arg, real_t* res) {
   *   <FUNCTION_BODY>
   * }
   *
   */
  CASADI_EXPORT Function jit(const std::string& name, int n_in, int n_out,
                             const std::string& body, const Dict& opts=Dict());

#ifndef SWIG
  /** \brief Create a just-in-time compiled function from a .casadi file
   */
  CASADI_EXPORT Function jit(const ParsedFile& file);
#endif // SWIG

} // namespace casadi

#endif // CASADI_JIT_HPP
