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

#ifndef SX_C_HPP
#define SX_C_HPP

#include "sx_c.h"
#include "../sx/sx.hpp"
#include "c_interface.hpp"

namespace CasADi{

/// Convert a pointer to an SX to a reference to an SX
SX& sx_ref(sx_ptr ptr);

/// Convert a pointer to an vector<SX> to a reference to an vector<SX>
std::vector<SX>& sx_vec(sx_vec_ptr v);

}

#endif // SX_C_HPP
