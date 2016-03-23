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


#include "casadi/core/polynomial.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/function/function.hpp"
#include "direct_collocation.hpp"

using namespace std;
namespace casadi {

  DirectCollocation::DirectCollocation(const Dict& f)
      : Function(  ) {
    addOption("interpolation_order",           OT_INT,  3,
              "Order of the interpolating polynomials");
    addOption("collocation_scheme",            OT_STRING,  "radau",
              "Collocation scheme", "radau|legendre");
    addOption("map_strategy",         OT_STRING,  "unrolled",
              "Collocation scheme", "unrolled|symbolic_serial|symbolic_openmp");
    setOption("name", "unnamed_direct_collocation");
  }

  DirectCollocation::getNlp() {
  }

} // namespace casadi
