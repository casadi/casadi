/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef SX_MATRIX_C_HPP
#define SX_MATRIX_C_HPP

#include "sx_matrix_c.h"
#include "../sx/sx_matrix.hpp"
#include "c_interface.hpp"
#include "sx_c.hpp"

namespace CasADi{

  SXMatrix& get_sx_matrix(sx_matrix_ref ref);

  std::vector<SXMatrix>& get_sx_matrix_vec(sx_matrix_vec v);

}
  
#endif // SX_MATRIX_C_HPP
