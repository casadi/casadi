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

#include "matrix_element.hpp"
using namespace std;

namespace CasADi{

MatrixElement::MatrixElement(const MX& x, int i, int j) : i_(i), j_(j){
  setDependencies(x);
  setSparsity(CRSSparsity(1,1,true));
  
  // Only one non-zero
  int ind = x->sparsity().getNZ(i,j);
  casadi_assert(ind>=0);
  
  nzind_.resize(1,ind);
  argind_.resize(1,0);
}

MatrixElement* MatrixElement::clone() const{
  return new MatrixElement(*this);
}

void MatrixElement::print(std::ostream &stream) const{
  stream << dep(0) << "(" << i_ << "," << j_ << ")";
}

} // namespace CasADi
