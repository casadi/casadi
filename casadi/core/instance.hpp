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


#ifndef CASADI_INSTANCE_HPP
#define CASADI_INSTANCE_HPP

#include "layout.hpp"
#include <vector>

namespace casadi {

  // Move to internal
  struct Instance {
    std::vector<bool> arg_null;
    std::vector<bool> res_null;
    std::vector<casadi_int> stride_in;
    std::vector<casadi_int> stride_out;
    bool prefer_inline;
    bool operator==(const Instance &rhs) const { return rhs.arg_null==arg_null && rhs.res_null==res_null && stride_in ==rhs.stride_in && stride_out==rhs.stride_out && prefer_inline==rhs.prefer_inline; }
  };

} // namespace casadi

#endif // CASADI_INSTANCE_HPP
