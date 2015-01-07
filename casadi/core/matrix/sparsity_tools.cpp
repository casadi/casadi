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


#include "sparsity_tools.hpp"
#include "../casadi_exception.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace casadi {

  int rank(const Sparsity& a) {
    std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    a.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
    return coarse_colblock.at(3);
  }

} // namespace casadi
