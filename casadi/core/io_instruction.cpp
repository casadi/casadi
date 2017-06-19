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


#include "io_instruction.hpp"

using namespace std;

namespace casadi {
  Input::Input(const Sparsity& sp, int ind, int segment, int offset)
    : IOInstruction(ind, segment, offset) {
    set_sparsity(sp);
  }

  std::string Input::print(const std::vector<std::string>& arg) const {
    stringstream s;
    s << "input[" << ind_ << "][" << segment_ << "]";
    return s.str();
  }

  Output::Output(const MX& x, int ind, int segment, int offset)
    : IOInstruction(ind, segment, offset) {
    set_dep(x);
  }

  std::string Output::print(const std::vector<std::string>& arg) const {
    stringstream s;
    s << "output[" << ind_ << "][" << segment_ << "]";
    return s.str();
  }

} // namespace casadi
